"""Expose the existing exp-transformed Laplace draw in a separate scalar logger."""
import re


OLD_SIGNATURE = 'model sequenceData referenceData initializationAudit tree logParamsTSV logParamsJSON [logA] [logCatStates] ='
NEW_SIGNATURE = OLD_SIGNATURE.replace('logParamsJSON [logA]', 'logParamsJSON logAlphaStateV8 [logA]')
OLD_CALL = '(model sequenceData referenceData isTest tree logParamsTSV logParamsJSON [logA] [logCatStates])'
NEW_CALL = OLD_CALL.replace('logParamsJSON [logA]', 'logParamsJSON logAlphaStateV8 [logA]')
OLD_BINDING = ';(smodel,log_smodel) <- sample_smodel aa\n'
NEW_BINDING = ';(smodel,log_smodel,projectAlphaFieldsV8) <- sample_smodel aa\n'
ALPHA_SETUP = ';logAlphaStateV8 <- if jsonEnabled then projectScalarJSONLoggerV6 (outputDirectory </> "C1.P1.log-alpha-samples.jsonl") else return noLogger\n'
ALPHA_ACTION = ';addLogger $ (logAlphaStateV8 (LoggerValues ["S1" %>% projectAlphaFieldsV8] (contextLogValues loggerValues)))\n'
DRAW = re.compile(r';alpha_2 <- sample \(logLaplace ([0-9]+) ([0-9]+)\)\n')
LATENT_DRAW = re.compile(r';projectLatentLogAlphaV8 <- sample \(laplace ([0-9]+) ([0-9]+)\)\n;let \{alpha_2 = exp projectLatentLogAlphaV8\}\n')


def fields(mu, scale):
    return (';let {projectAlphaFieldsV8 = ["latentLogAlpha" %=% projectLatentLogAlphaV8,'
            '"derivedAlpha" %=% alpha_2,"categoryRates" %=% (values (SModel.gammaRatesMean alpha_2 4)),'
            '"laplaceLocation" %=% (' + mu + ' :: Double),"laplaceScale" %=% (' + scale + ' :: Double),'
            '"latentLogDensity" %=% (ln (pdf (laplace ' + mu + ' ' + scale + ') projectLatentLogAlphaV8))]}\n')


def transform(text):
    assert 'projectLatentLogAlphaV8' not in text and '-- BEGIN project joint FASTA lines v7' in text
    matches = list(DRAW.finditer(text))
    assert len(matches) == 1
    mu, scale = matches[0].groups()
    assert (mu, scale) in [('0', '2'), ('0', '1'), ('6', '2')]
    start = text.index('sample_smodel alpha  =')
    end = text.index('\nsample_imodel ', start)
    section = text[start:end]
    assert section.count(';return (result, loggers)\n') == 1
    section = DRAW.sub(';projectLatentLogAlphaV8 <- sample (laplace ' + mu + ' ' + scale + ')\n'
                       ';let {alpha_2 = exp projectLatentLogAlphaV8}\n', section)
    section = section.replace(';return (result, loggers)\n',
                              fields(mu, scale) + ';return (result, loggers, projectAlphaFieldsV8)\n')
    text = text[:start] + section + text[end:]
    for before, after in [(OLD_SIGNATURE, NEW_SIGNATURE), (OLD_CALL, NEW_CALL), (OLD_BINDING, NEW_BINDING)]:
        assert text.count(before) == 1
        text = text.replace(before, after)
    anchor = ';logA <- if loggingEnabled then alignmentLogger '
    assert text.count(anchor) == 1
    text = text.replace(anchor, ALPHA_SETUP + anchor)
    anchor = ';addLogger $ (logParamsJSON loggerValues)\n'
    assert text.count(anchor) == 1
    return text.replace(anchor, anchor + ALPHA_ACTION)


def reverse(text):
    matches = list(LATENT_DRAW.finditer(text))
    assert len(matches) == 1
    mu, scale = matches[0].groups()
    text = LATENT_DRAW.sub(';alpha_2 <- sample (logLaplace ' + mu + ' ' + scale + ')\n', text)
    marker = fields(mu, scale) + ';return (result, loggers, projectAlphaFieldsV8)\n'
    assert text.count(marker) == 1
    text = text.replace(marker, ';return (result, loggers)\n')
    for before, after in [(NEW_SIGNATURE, OLD_SIGNATURE), (NEW_CALL, OLD_CALL), (NEW_BINDING, OLD_BINDING),
                          (ALPHA_SETUP, ''), (ALPHA_ACTION, '')]:
        assert text.count(before) == 1
        text = text.replace(before, after)
    assert 'projectLatentLogAlphaV8' not in text and 'logAlphaStateV8' not in text
    return text
