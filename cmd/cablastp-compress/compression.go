package main

import (
  "bytes"
  "runtime"
  "sync"

  "github.com/ndaniels/cablastp2"
)

// compressPool represents a pool of workers where each worker is responsible
// for compressing a single sequence at a time.
type compressPool struct {
  db     *cablastp.DB
  jobs   chan compressJob
  wg     *sync.WaitGroup
  closed bool
}

// compressJob values are messages sent to the pool of workers when a new
// sequence should be compressed.
type compressJob struct {
  orgSeqId int
  orgSeq   *cablastp.OriginalSeq
}

// startCompressWorkers initializes a pool of compression workers.
//
// The compressPool returned can be used to compress sequences concurrently.
func startCompressWorkers(db *cablastp.DB) compressPool {
  wg := &sync.WaitGroup{}
  jobs := make(chan compressJob, 200)
  pool := compressPool{
    db:     db,
    jobs:   jobs,
    wg:     wg,
    closed: false,
  }
  for i := 0; i < max(1, runtime.GOMAXPROCS(0)); i++ {
    wg.Add(1)
    go pool.worker()
  }
  return pool
}

// compress will construct a compressJob and send it to the worker pool.
//
// compress returns the next original sequence id to be used.
func (pool compressPool) compress(id int, seq *cablastp.OriginalSeq) int {
  pool.jobs <- compressJob{
    orgSeqId: id,
    orgSeq:   seq,
  }
  return id + 1
}

// worker is meant to be run as a goroutine. It allocates a goroutine-specific
// memory arena (to prevent allocation in hot spots like alignment and
// seed lookup), and sends the compressed sequences to the compressed
// database for writing.
func (pool compressPool) worker() {
  mem := newMemory()
  for job := range pool.jobs {
    comSeq := compress(pool.db, job.orgSeqId, job.orgSeq, mem)
    pool.db.ComDB.Write(comSeq)
  }
  pool.wg.Done()
}

// done 'joins' the worker goroutines. (Blocks until all workers are finished
// compressing sequences.)
func (pool *compressPool) done() {
  if pool.closed {
    return
  }
  pool.closed = true
  close(pool.jobs)
  pool.wg.Wait()
}

// TODO: need to pass around a pair of nativeSeq and reducedSeq
// What gets saved to compressed db is nativeSeq
// What gets saved to coarse fasta is reducedSeq

// compress will convert an original sequence into a compressed sequence.
// The process involves finding commonality in the original sequence with
// other sequences in the coarse database, and linking those common
// sub-sequences to sub-sequences in the coarse database.
//
// N.B. `mem` is used in alignment and seed lookups to prevent allocation.
// Think of them as goroutine-specific memory arenas.
func compress(db *cablastp.DB, orgSeqId int,
  orgSeq *cablastp.OriginalSeq, mem *memory) cablastp.CompressedSeq {

  // cseqExt and oseqExt will contain `extSeedSize` residues after the end
  // of any particular seed in coarse and original sequences, respectively.
  // If the residues are not equivalent, that particular seed is skipped.
  var cseqExt, oseqExt []byte

  // redSeq will hold the reduced-alphabet sequence
  redSeq := cablastp.NewReducedSeq(orgSeq)

  // Start the creation of a compressed sequence.
  cseq := cablastp.NewCompressedSeq(orgSeqId, orgSeq.Name)

  // Convenient aliases.
  coarsedb := db.CoarseDB
  mapSeedSize := db.MapSeedSize
  extSeedSize := db.ExtSeedSize
  olen := orgSeq.Len()

  // Keep track of two pointers. 'current' refers to the residue index in the
  // original sequence that extension is currently originating from.
  // 'lastMatch' refers to the residue index of the *end* of the last match
  // with a coarse sequence in the compressed database.
  lastMatch, current := 0, 0

  // Iterate through the original sequence a 'kmer' at a time.
  skipSize := 4
  limit := olen-mapSeedSize-extSeedSize-skipSize
  for current = 0; current <= limit; current += skipSize {
    kmer := redSeq.Residues[current : current+mapSeedSize]

    // skip wildcard-containing kmers
    if bytes.IndexByte(kmer, 'N') > -1 {
      continue
    }

    seeds := coarsedb.Seeds.Lookup(kmer, &mem.seeds)

    // Before trying to extend this with seeds, check to see if there is
    // a low complexity region within `db.MinMatchLen` residues from
    // `current`. If there is, skip ahead to the end of it.
    // note that we are checking the original, unreduced sequence here!
    if db.LowComplexity > 0 {
      skip := skipLowComplexity(
        orgSeq.Residues[current:], db.MinMatchLen, db.LowComplexity)
      if skip > 0 {
        current += skip
        continue
      }
    }

    // Each seed location corresponding to the current K-mer must be
    // used to attempt to extend a match.
    for _, seedLoc := range seeds {
      corSeqId := int(seedLoc[0])
      corResInd := int(seedLoc[1])
      corSeq := coarsedb.CoarseSeqGet(uint(corSeqId))

      // If the seed extension extends beyond the end of the coarse
      // sequence pointed to by seedLoc, then move along.
      extCorStart := corResInd + mapSeedSize
      extOrgStart := current + mapSeedSize
      if extCorStart+extSeedSize >= corSeq.Len() {
        continue
      }

      // If the seed extensions in each sequence are not equivalent,
      // skip this seedLoc.
      cseqExt = corSeq.Residues[extCorStart : extCorStart+extSeedSize]
      oseqExt = redSeq.Residues[extOrgStart : extOrgStart+extSeedSize]
      if !bytes.Equal(cseqExt, oseqExt) {
        continue
      }

      // The "match" between coarse and original (reduced) sequence will
      // occur somewhere between the the residue index of the seed and
      // the end of the sequence for the coarse sequence, and the
      // position of the "current" pointer and the end of the sequence
      // for the original sequence.
      // TODO maybe the simplest way to match backwards is to reverse the
      // coarse and reduced sequences FROM the last match TO the current,
      // and call extendMatch on those as well. Reverse the result, and
      // prepend it to corMatch and redMatch here.
      corMatch, redMatch := extendMatch(
        corSeq.Residues[corResInd:], redSeq.Residues[current:],
        db.GappedWindowSize, db.UngappedWindowSize,
        db.MatchKmerSize, db.ExtSeqIdThreshold,
        mem)
        
      // potentially extend this match back as far as the lastMatch (for redSeq)
      // and beginning of the corSeq
      backCorMatch, backRedMatch := extendMatch(
        reverse(corSeq.Residues[0:corResInd]), 
        reverse(redSeq.Residues[lastMatch:current]),
        db.GappedWindowSize, db.UngappedWindowSize,
        db.MatchKmerSize, db.ExtSeqIdThreshold,
        mem)
        
      redMatch = append(reverse(backRedMatch), redMatch...)
      corMatch = append(reverse(backCorMatch), corMatch...)
      current -= len(redMatch)
      corResInd -= len(corMatch)

      // If the part of the original (reduced) sequence does not exceed the
      // minimum match length, then we don't accept the match and move
      // on to the next one.
      if len(redMatch) < db.MinMatchLen {
        continue
      }

      // TODO: we don't need to realign here if we are just
      // preserving original sequences.

      // alignment := nwAlign(corMatch, orgMatch, mem)
      // 	id := cablastp.SeqIdentity(alignment[0], alignment[1])
      // 	if id < db.MatchSeqIdThreshold {
      // 		continue
      // 	}

      // If we end up extending a match because we're close to
      // some boundary (either a sequence or a match boundary), then
      // we need to perform another alignment.
      // changed := false

      // If we're close to the end of the original sequence, extend
      // the match to the end.
      if len(redMatch)+db.MatchExtend >= redSeq.Len()-int(current) {
        redMatch = redSeq.Residues[current:]
        // changed = true
      }

      // And if we're close to the end of the last match, extend this
      // match backwards.
      if current-lastMatch <= db.MatchExtend {
        end := current + len(redMatch)
        redMatch = redSeq.Residues[lastMatch:end]
        current = lastMatch
        // changed = true
      }

      // Otherwise, we accept the first valid match and move on to the
      // next kmer after the match ends.
      corStart := corResInd
      corEnd := corStart + len(corMatch)
      orgStart := current
      orgEnd := orgStart + len(redMatch)

      // If there are residues between the end of the last match
      // and the start of this match, then that means no good match
      // could be found for those residues. Thus, they are added to
      // the coarse database. (A pathological LinkToCoarse is
      // created with an empty diff script that points to the added
      // region in the coarse database in its entirety.)
      if orgStart-lastMatch > 0 {
        redSub := redSeq.NewSubSequence(
          uint(lastMatch), uint(current))
        orgSub := orgSeq.NewSubSequence(
          uint(lastMatch), uint(current))
        addWithoutMatch(&cseq, coarsedb, orgSeqId, orgSub, redSub)
      }

      // For the given match, add a LinkToCoarse to the portion of
      // the coarse sequence matched. This serves as a component
      // of a compressed original sequence. Also, add a
      // LinkToCompressed to the coarse sequence matched. This
      // serves as a bridge to expand coarse sequences into their
      // original sequences.
      orgMatch := string(orgSeq.Residues[orgStart:orgEnd])
      cseq.Add(cablastp.NewLinkToCoarse(
        uint(corSeqId), uint(corStart), uint(corEnd), orgMatch))
      corSeq.AddLink(cablastp.NewLinkToCompressed(
        uint32(orgSeqId), uint16(corStart), uint16(corEnd)))

      // Skip the current pointer ahead to the end of this match.
      // Update the lastMatch pointer to point at the end of this
      // match.
      lastMatch = orgEnd
      current = orgEnd - 1

      // Don't process any more seedLocs for this K-mer once we've
      // found a match.
      break
    }
  }

  // If there are any leftover residues, then no good match for them
  // could be found. Therefore, add them to the coarse database and
  // create the appropriate links.
  if orgSeq.Len()-lastMatch > 0 {
    orgSub := orgSeq.NewSubSequence(uint(lastMatch), uint(orgSeq.Len()))
    redSub := redSeq.NewSubSequence(uint(lastMatch), uint(redSeq.Len()))
    addWithoutMatch(&cseq, coarsedb, orgSeqId, orgSub, redSub)
  }

  return cseq
}

func reverse(a []byte) []byte {
  l := len(a)
  result := make([]byte, l)
  for i, v := range a {
    result[l-i-1] = v
  }
  return result
}

// extendMatch uses a combination of ungapped and gapped extension to find
// quality candidates for compression.
func extendMatch(corRes, orgRes []byte,
  gappedWindowSize, ungappedWindowSize, kmerSize, idThreshold int,
  mem *memory) (corMatchRes, orgMatchRes []byte) {

  // Starting at seedLoc.resInd and current (from 'compress'), corMatchLen
  // and orgMatchLen correspond to the length of the match in each of
  // the coarse and the original sequence, respectively.
  // At the end of the loop, the slices [seedLoc.resInd:corMatchLen]
  // and [current:orgMatchLen] will correspond to the match. (Again, this
  // is in the context of the inner loop in 'compress'. For this particular
  // function, corMatchLen and orgMatch start at 0, so that the matches
  // eventually returned correspond to the [:corMatchLen] and [:orgMatchLen]
  // slices.)
  corMatchLen, orgMatchLen := 0, 0
  for {
    // If the match has consumed either of the coarse or original
    // sequence, then we must quit with what we have.
    if corMatchLen == len(corRes) || orgMatchLen == len(orgRes) {
      break
    }

    // Ungapped extension returns an integer corresponding to the
    // number of residues that the match was extended by.
    matchLen := alignUngapped(
      corRes[corMatchLen:], orgRes[orgMatchLen:],
      ungappedWindowSize, kmerSize, idThreshold)

    // Since ungapped extension increases the coarse and
    // original sequence match portions equivalently, add the
    // match length to both.
    corMatchLen += matchLen
    orgMatchLen += matchLen

    // Gapped extension returns an alignment corresponding to the
    // window starting after the previous ungapped extension
    // ended plus the gapped window size. (It is bounded by the
    // length of each sequence.)
    alignment := nwAlign(
      corRes[corMatchLen:min(len(corRes), corMatchLen+gappedWindowSize)],
      orgRes[orgMatchLen:min(len(orgRes), orgMatchLen+gappedWindowSize)],
      mem)

    // If the alignment has a sequence identity below the
    // threshold, then gapped extension has failed. We therefore
    // quit and are forced to be satisfied with whatever
    // corMatchLen and orgMatchLen are set to.
    id := cablastp.SeqIdentity(alignment[0], alignment[1])
    if id < idThreshold {
      break
    }

    // We live to die another day.
    // We need to add to the corMatch{Pos,Len} and orgMatch{Pos,Len}
    // just like we did for ungapped extension. However, an
    // alignment can correspond to two different sized subsequences
    // of the coarse and original sequence. Therefore, only
    // increase each by the corresponding sizes from the
    // alignment.
    corMatchLen += alignLen(alignment[0])
    orgMatchLen += alignLen(alignment[1])
  }

  return corRes[:corMatchLen], orgRes[:orgMatchLen]
}

// skipLowComplexity looks for a low complexity region starting at the
// beginning of `seq` and up to `windowSize`. If one is found, `x` is returned
// where `x` corresponds to the position of the first residue after
// the low complexity region has ended. If a low complexity region isn't
// found, `0` is returned.
//
// N.B. regionSize is the number of contiguous positions in the sequence
// that must contain the same residue in order to qualify as a low complexity
// region.
func skipLowComplexity(seq []byte, windowSize, regionSize int) int {
  upto := min(len(seq), windowSize+regionSize)
  last, repeats, i, found := byte(0), 1, 0, false
  for i = 0; i < upto; i++ {
    if seq[i] == last {
      repeats++
      if repeats >= regionSize {
        found = true
        break
      }
      continue
    }

    // The last residue isn't the same as this residue, so reset.
    last = seq[i]
    repeats = 1
  }
  if !found { // no low complexity region was found.
    return 0
  }

  // We're in a low complexity region. Consume as many residues equal
  // to `last` as possible.
  //
  // N.B. `i` is already set to where we left off in the last loop.
  for ; i < len(seq); i++ {
    if seq[i] != last { // end of low complexity region
      break
    }
  }
  return i
}

// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the coarse database to the coarse database.
// A LinkToCompressed is created and automatically added to the new coarse
// sequence.
//
// An appropriate link is also added to the given compressed sequence.
// TODO we need to juggle the reduced and original seq portions
func addWithoutMatch(cseq *cablastp.CompressedSeq,
  coarsedb *cablastp.CoarseDB, orgSeqId int, orgSub *cablastp.OriginalSeq,
  redSub *cablastp.ReducedSeq) {

  // Explicitly copy residues to avoid pinning memory.
  redSubCpy := make([]byte, len(redSub.Residues))
  copy(redSubCpy, redSub.Residues)

  corSeqId, corSeq := coarsedb.Add(redSubCpy)
  corSeq.AddLink(
    cablastp.NewLinkToCompressed(uint32(orgSeqId), 0, uint16(len(redSubCpy))))

  cseq.Add(
    cablastp.NewLinkToCoarse(uint(corSeqId), 0, uint(len(redSubCpy)),
      string(orgSub.Residues)))
}

func min(a, b int) int {
  if a < b {
    return a
  }
  return b
}

func max(a, b int) int {
  if a > b {
    return a
  }
  return b
}
