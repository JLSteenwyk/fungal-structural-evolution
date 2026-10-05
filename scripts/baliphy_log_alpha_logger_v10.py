"""Retain the installed exp-transform's nested Random fmap while capturing its draw."""
import re

import baliphy_log_alpha_logger_v9 as v9


DIRECT = re.compile(r';projectLatentLogAlphaV10 <- sample \(laplace ([0-9]+) ([0-9]+)\)\n;let \{alpha_2 = exp projectLatentLogAlphaV10\}\n')
CAPTURE = re.compile(r';\(alpha_2,projectLatentLogAlphaV10\) <- \(\\projectDrawV10 -> \(exp projectDrawV10,projectDrawV10\)\) <\$> sample \(laplace ([0-9]+) ([0-9]+)\)\n')


def fields(mu,scale):
    return v9.fields(mu,scale).replace('V9','V10')


def transform(text):
    result = v9.transform(text).replace('V9','V10')
    matches = list(DIRECT.finditer(result))
    assert len(matches) == 1
    mu,scale = matches[0].groups()
    replacement = (';(alpha_2,projectLatentLogAlphaV10) <- '
                   '(\\projectDrawV10 -> (exp projectDrawV10,projectDrawV10)) <$> sample '
                   '(laplace ' + mu + ' ' + scale + ')\n')
    return result[:matches[0].start()] + replacement + result[matches[0].end():]


def reverse(text):
    matches = list(CAPTURE.finditer(text))
    assert len(matches) == 1
    mu,scale = matches[0].groups()
    direct = (';projectLatentLogAlphaV10 <- sample (laplace ' + mu + ' ' + scale + ')\n'
              ';let {alpha_2 = exp projectLatentLogAlphaV10}\n')
    result = text[:matches[0].start()] + direct + text[matches[0].end():]
    return v9.reverse(result.replace('V10','V9'))
