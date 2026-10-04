"""Reversibly wrap logger IO actions; leave model expressions and logger bodies intact."""

IMPORT = 'import MCMC.Loggers (LoggerAction)\n'
BLOCK = '''-- BEGIN project diagnostic logger phase trace v1
projectTracePhase :: String -> String -> String -> Int -> IO ()
projectTracePhase name phase boundary iteration = do
  hPutStrLn stderr ("PROJECT_LOGGER_PHASE " ++ name ++ " " ++ phase ++ " " ++ boundary ++ " " ++ show iteration)
  hFlush stderr

projectTraceLogger :: String -> LoggerAction -> LoggerAction
projectTraceLogger name (getFields, writeSample) = (tracedFields, tracedWrite)
  where
    tracedFields iteration context = do
      projectTracePhase name "context" "begin" iteration
      fields <- getFields iteration context
      projectTracePhase name "context" "end" iteration
      return fields
    tracedWrite iteration fields = do
      projectTracePhase name "write" "begin" iteration
      writeSample iteration fields
      projectTracePhase name "write" "end" iteration
-- END project diagnostic logger phase trace v1

'''
LABELS = ['scalar-tsv', 'scalar-cjson', 'ancestral-alignment', 'joint-states']


def transform(text):
    assert IMPORT not in text and BLOCK not in text
    assert text.count('import Probability.Logger\n') == 1
    text = text.replace('import Probability.Logger\n', 'import Probability.Logger\n' + IMPORT)
    anchor = '-- BEGIN project scalar JSON logger v6\n'
    assert text.count(anchor) == 1
    text = text.replace(anchor, BLOCK + anchor)
    rows = text.splitlines(keepends=True)
    indices = [i for i, row in enumerate(rows) if row.startswith(';addLogger $ ')]
    assert len(indices) == len(LABELS)
    for i, name in zip(indices, LABELS):
        rows[i] = rows[i].replace(';addLogger $ ', ';addLogger $ projectTraceLogger "' + name + '" $ ', 1)
    return ''.join(rows)


def reverse(text):
    assert text.count(IMPORT) == text.count(BLOCK) == 1
    text = text.replace(IMPORT, '').replace(BLOCK, '')
    for name in LABELS:
        traced = ';addLogger $ projectTraceLogger "' + name + '" $ '
        assert text.count(traced) == 1
        text = text.replace(traced, ';addLogger $ ')
    return text


def events(text):
    result = []
    for line in text.splitlines():
        if not line.startswith('PROJECT_LOGGER_PHASE '):
            continue
        tag, name, phase, boundary, iteration = line.split()
        assert name in LABELS and phase in ['context', 'write'] and boundary in ['begin', 'end']
        result.append(dict(logger=name, phase=phase, boundary=boundary, iteration=int(iteration)))
    return result
