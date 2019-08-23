import { dnaComplement } from "../../Utils/parser";
import {
  calcTm,
  getMismatchIndices,
  returnRanges,
  reverse
} from "../../Utils/sequence";

/**
 * Given a list of primer mismatches and the primer's annealing sequence
 * Returns the mismatch indices and the part of primer sequence that actually anneals
 * @param {string} sequence
 * @param {string} subVector
 * @return {Array} mismatches as an array of array of start/end indices
 * @return {string} the string sequence of the part of the primer that anneals
 */
const findMismatches = (sequence, subVector) => {
  let annealSequence = sequence;
  let mismatches = returnRanges(getMismatchIndices(sequence, subVector));

  // Logic for when to combine mismatches into an overhang tail
  // If the remaining sequence after a mismatch is more than 25% mismatch a tail is formed
  if (mismatches.length > 0) {
    let index = 0;
    while (index < mismatches.length) {
      const remainingMismatches = mismatches.slice(
        0,
        mismatches.length - index
      );
      if (remainingMismatches.length < 2) {
        break;
      }
      const mismatchTotalLength = remainingMismatches.reduce(
        (acc, mismatch) => acc + (mismatch[1] + 1 - mismatch[0]),
        0
      );
      const mismatchAreaLength = sequence.slice(
        0,
        remainingMismatches[remainingMismatches.length - 1][1] + 1
      ).length;

      if (mismatchTotalLength / mismatchAreaLength > 0.25) {
        mismatches = mismatches
          .slice(mismatches.length - index, mismatches.length)
          .concat([[0, mismatches[mismatches.length - 1 - index][1]]]);

        annealSequence = sequence.slice(
          mismatches[mismatches.length - 1][1] + 1
        );
        break;
      }
      index += 1;
    }
    annealSequence = sequence.slice(mismatches[mismatches.length - 1][1] + 1);
  }
  mismatches = mismatches.map(mismatch => ({
    start: mismatch[0],
    end: mismatch[1] + 1 // because mismatches return indices of mismatch and we want to bound to end after the last index
  }));

  return { mismatches, annealSequence };
};

/**
 * Find binding sites on one strand of DNA
 * @param {Array} primers
 * @param {string} vectorSequence
 * @param {string} direction
 * @return {Array} array of primers with viewer meta information
 */
const findBindingSites = (primers = [], vectorSeq, direction) => {
  const matchLength = 10;
  const minTm = 40;

  const primerBindingSites = [];
  const forward = direction === "FORWARD";

  primers.forEach((primer, i) => {
    const { overhang = "" } = primer;
    let { sequence, strict } = primer;
    strict = strict || false;
    if (sequence === "") {
      return;
    }
    sequence = sequence.toLowerCase();
    const sequenceLength = sequence.length;
    console.log(`BINDING CHECK PRIMERSEQ${i}: `, sequence);
    console.log(
      `BINDING ${
        sequenceLength >= 10 ? "CHECK" : "FAIL"
      } PRIMERSEQ${i} length ${sequenceLength} greater than or equal to 10`
    );
    const vectorSequence = vectorSeq.toLowerCase();
    const vectorLength = vectorSequence.length;

    const expandedVectorSequence =
      vectorSequence + vectorSequence.substring(0, sequenceLength); // Used for looking for binding sites that cross 0 index

    let annealSequence = sequence;
    let { mismatches, matchSeq } = [];

    matchSeq =
      sequenceLength < matchLength
        ? sequence
        : sequence.substring(sequenceLength - matchLength, sequenceLength);

    matchSeq = forward ? matchSeq : reverse(matchSeq);
    console.log(
      `BINDING ${
        matchSeq.length >= 10 ? "CHECK" : "FAIL"
      } PRIMERMATCHSEQ${i} ${matchSeq} length greater than or equal to 10`
    );

    const regex = new RegExp(matchSeq, "gi");
    let result = regex.exec(expandedVectorSequence);
    const combinedSequence = (overhang || "").concat(sequence);
    while (result) {
      if (result.index < vectorLength) {
        const tailCrossZero = forward
          ? result.index + matchSeq.length - sequenceLength < 0
          : result.index + sequenceLength > vectorLength;
        const headCrossZero = forward
          ? vectorLength - (result.index + matchSeq.length) < 0
          : result.index + matchSeq.length > vectorLength;
        const crossZero = tailCrossZero || headCrossZero;

        let startIndex = forward
          ? result.index - sequenceLength + matchSeq.length
          : result.index;
        let endIndex = startIndex + sequenceLength;
        let subVector = vectorSequence.substring(startIndex, endIndex);

        if (crossZero) {
          if (forward) {
            startIndex = tailCrossZero
              ? vectorLength - (sequenceLength - result.index - matchSeq.length)
              : result.index - sequenceLength + matchSeq.length;
          }
          endIndex = sequenceLength - (vectorLength - startIndex);
          subVector =
            vectorSequence.substring(startIndex, vectorLength) +
            vectorSequence.substring(0, endIndex);
        }

        subVector = forward ? subVector : reverse(subVector);

        const matchTm = calcTm(sequence, subVector);
        console.log(
          `BINDING ${
            matchTm >= primer.tm || matchTm > minTm ? "CHECK" : "FAIL"
          } SUBVECTOR${i} ${subVector} has Tm ${matchTm} greater than or equal to ${
            primer.tm
          } or greater than ${minTm}`
        );
        console.log(
          `BINDING ${
            (forward && endIndex < vectorLength) || startIndex < vectorLength
              ? "CHECK"
              : "FAIL"
          } SUBVECTOR${i} has correct directionality`
        );
        console.log(
          `BINDING ${
            result.index < vectorLength ? "CHECK" : "FAIL"
          } SUBVECTOR${i} has index within vector length`
        );
        if (
          (matchTm >= primer.tm || matchTm > minTm) &&
          ((forward && endIndex < vectorLength) || startIndex < vectorLength) &&
          result.index < vectorLength
        ) {
          if (overhang) {
            if (forward) {
              startIndex -= overhang.length;
              if (startIndex < 0) startIndex = vectorLength + startIndex;
            }
            if (!forward) {
              endIndex += overhang.length;
              if (endIndex > vectorLength) endIndex -= vectorLength;
            }
            ({ mismatches = [], annealSequence = "" } = findMismatches(
              combinedSequence,
              "X".repeat(overhang.length).concat(subVector)
            ));
            if (mismatches[0] && mismatches[0].start - overhang.length === 0) {
              mismatches[0].start = 0;
            } else {
              mismatches.push({ start: 0, end: overhang.length });
            }
          } else {
            ({ mismatches = [], annealSequence = "" } = findMismatches(
              sequence,
              subVector
            ));
          }
          console.log(
            `BINDING CHECK MISMATCHES${i} ${mismatches}, ANNEALING${i} ${annealSequence}`
          );
          const uniqMismatch = {};
          mismatches = mismatches
            .sort((a, b) => a.start < b.start)
            .map(m => ({
              id: `${m.start}-${m.end}`,
              start: m.start,
              end: m.end
            }))
            .filter(m => {
              if (uniqMismatch[m.id]) {
                return false;
              }
              uniqMismatch[m.id] = true;
              return true;
            })
            .map(m => ({ start: m.start, end: m.end }));
          console.log(
            `BINDING CHECK ${strict} should be false, BINDING COUNT AT ITERATION ${i} is ${
              primerBindingSites.length
            }`
          );
          if ((strict && mismatches.length < 1) || !strict) {
            primerBindingSites.push({
              ...primer,
              id: `${primer.id}-${startIndex}`,
              sequence: combinedSequence,
              start: startIndex,
              end: endIndex,
              direction: direction,
              mismatches: mismatches,
              annealSequence: annealSequence,
              strict: strict
            });
          }
        }
      }
      result = regex.exec(expandedVectorSequence);
    }
  });
  return primerBindingSites;
};

/**
 * Gives primers meta information needed by sequence viewers
 * @param {Array} primers: array of primers
 * @param {string} vector
 * @return {Array} array of primers with viewer meta information
 */
export const findAllBindingSites = (primers, vector) => {
  const { seq: vectorSeq, compSeq: vectorComp } = dnaComplement(vector);

  const bindings = findBindingSites(primers, vectorSeq, "FORWARD").concat(
    findBindingSites(primers, vectorComp, "REVERSE")
  );
  console.log("BINDINGS FINAL: ", primers, bindings);
  return bindings;
};

export default findAllBindingSites;
