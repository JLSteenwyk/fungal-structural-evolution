"""Construct joint FASTA lines from native sequence rows, avoiding whole-alignment Char lists."""

BLOCK = '''-- BEGIN project joint FASTA lines v7
projectJointFastaLines characterData = concat
  [[Text.cons '>' label,
    sequenceToText (getAlphabet characterData) (getAmbiguities characterData) sequence]
   | (label, sequence) <- getSequences characterData]
-- END project joint FASTA lines v7

'''
OLD_ALIGNMENT = ';let {ancAlignment  = toFasta (ancestralAlignment tree alignment (getSMap smodel) aa ancStates)}\n'
NEW_ALIGNMENT = (';let {projectAncestralDataV7 = ancestralAlignment tree alignment (getSMap smodel) aa ancStates}\n'
                 ';let {ancAlignment  = toFasta projectAncestralDataV7}\n')
OLD_LINES = '(map Text.pack (lines (Text.unpack ancAlignment)))'
NEW_LINES = '(projectJointFastaLines projectAncestralDataV7)'


def transform(text):
    assert BLOCK not in text and text.count(OLD_ALIGNMENT) == text.count(OLD_LINES) == 1
    anchor = '-- BEGIN project scalar JSON logger v6\n'
    assert text.count(anchor) == 1
    return text.replace(anchor, BLOCK + anchor).replace(OLD_ALIGNMENT, NEW_ALIGNMENT).replace(OLD_LINES, NEW_LINES)


def reverse(text):
    assert text.count(BLOCK) == text.count(NEW_ALIGNMENT) == text.count(NEW_LINES) == 1
    return text.replace(BLOCK, '').replace(NEW_ALIGNMENT, OLD_ALIGNMENT).replace(NEW_LINES, OLD_LINES)
