"""Capture the same latent draw through the existing scalar logging action."""
import re
import json

from baliphy_log_alpha_logger_v8 import fields as v8_fields


DRAW = re.compile(r';alpha_2 <- sample \(logLaplace ([0-9]+) ([0-9]+)\)\n')
LATENT_DRAW = re.compile(r';projectLatentLogAlphaV9 <- sample \(laplace ([0-9]+) ([0-9]+)\)\n;let \{alpha_2 = exp projectLatentLogAlphaV9\}\n')
OLD_BINDING = ';(smodel,log_smodel) <- sample_smodel aa\n'
NEW_BINDING = ';(smodel,log_smodel,projectAlphaFieldsV9) <- sample_smodel aa\n'
OLD_SETUP = ';logParamsJSON <- if jsonEnabled then projectScalarJSONLoggerV6 (outputDirectory </> "C1.log.json") else return noLogger\n'
NEW_SETUP = ';logParamsJSON <- if jsonEnabled then projectLogAlphaLoggerV9 (outputDirectory </> "C1.log.json") (outputDirectory </> "C1.P1.log-alpha-samples.jsonl") else return noLogger\n'
OLD_ACTION = ';addLogger $ (logParamsJSON loggerValues)\n'
NEW_ACTION = ';addLogger $ (logParamsJSON (loggerValues,projectAlphaFieldsV9))\n'
BLOCK = '''-- BEGIN project latent log alpha logger v9
projectLogAlphaLoggerV9 filename alphaFilename = do
  handle <- openFile filename WriteMode
  alphaHandle <- openFile alphaFilename WriteMode
  let header = HEADER_PLACEHOLDER
  hPutStrLn handle header
  hPutStrLn alphaHandle header
  hFlush handle
  hFlush alphaHandle
  return (\\(loggerValues,alphaFields) ->
    (projectV6EvaluateContext (contextLogValues loggerValues),
     \\iteration contextFields -> do
       projectV6LogSample handle (parameterLogValues loggerValues) iteration contextFields
       projectV6LogSample alphaHandle ["S1" %>% alphaFields] iteration contextFields))
-- END project latent log alpha logger v9

'''
BLOCK = BLOCK.replace('HEADER_PLACEHOLDER', json.dumps('{"fields":["iter","prior","likelihood","posterior"],"nested":true,"format":"MCON","version":"0.2","projectScalarSchema":"native-cjson-explicit-special-values-v6"}'))


def fields(mu, scale):
    return v8_fields(mu, scale).replace('V8', 'V9')


def transform(text):
    assert 'projectLatentLogAlphaV9' not in text and '-- BEGIN project joint FASTA lines v7' in text
    matches = list(DRAW.finditer(text))
    assert len(matches) == 1
    mu, scale = matches[0].groups()
    assert (mu, scale) in [('0','2'), ('0','1'), ('6','2')]
    start = text.index('sample_smodel alpha  =')
    end = text.index('\nsample_imodel ', start)
    section = text[start:end]
    assert section.count(';return (result, loggers)\n') == 1
    section = DRAW.sub(';projectLatentLogAlphaV9 <- sample (laplace ' + mu + ' ' + scale + ')\n'
                       ';let {alpha_2 = exp projectLatentLogAlphaV9}\n', section)
    section = section.replace(';return (result, loggers)\n',
                              fields(mu,scale) + ';return (result, loggers, projectAlphaFieldsV9)\n')
    text = text[:start] + BLOCK + section + text[end:]
    for before, after in [(OLD_BINDING,NEW_BINDING), (OLD_SETUP,NEW_SETUP), (OLD_ACTION,NEW_ACTION)]:
        assert text.count(before) == 1
        text = text.replace(before,after)
    return text


def reverse(text):
    matches = list(LATENT_DRAW.finditer(text))
    assert len(matches) == 1
    mu, scale = matches[0].groups()
    text = LATENT_DRAW.sub(';alpha_2 <- sample (logLaplace ' + mu + ' ' + scale + ')\n',text)
    marker = fields(mu,scale) + ';return (result, loggers, projectAlphaFieldsV9)\n'
    assert text.count(marker) == text.count(BLOCK) == 1
    text = text.replace(marker,';return (result, loggers)\n').replace(BLOCK,'')
    for before, after in [(NEW_BINDING,OLD_BINDING), (NEW_SETUP,OLD_SETUP), (NEW_ACTION,OLD_ACTION)]:
        assert text.count(before) == 1
        text = text.replace(before,after)
    assert 'projectLatentLogAlphaV9' not in text and 'projectLogAlphaLoggerV9' not in text
    return text
