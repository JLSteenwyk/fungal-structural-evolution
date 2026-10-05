{-# LANGUAGE OverloadedStrings #-}
module Main where
import qualified Data.JSON as J
import qualified Data.Text as Text
import qualified Data.Text.IO as T
import Compiler.RealFloat (isNaN, isInfinite)
import Probability
import SModel.ASRV (gammaRatesMean)

tag :: Double -> String
tag x
  | isNaN x = "nan"
  | isInfinite x && x > 0 = "positive_infinity"
  | isInfinite x = "negative_infinity"
  | otherwise = "finite"

number x = if isNaN x || isInfinite x then J.String (Text.pack ("__log_alpha_diag_v2__:" ++ tag x)) else J.toJSON x

probe :: (String, Double, Double) -> Double -> IO ()
probe (label, mu, scale) logAlpha = do
  let alpha = exp logAlpha :: Double
      rates = values (gammaRatesMean alpha 4)
      basePrior = ln (pdf (laplace mu scale) logAlpha)
      derivedDensity = ln (pdf (logLaplace mu scale) alpha)
      row = J.Object [("prior", J.toJSON label), ("mu", J.toJSON mu), ("scale", J.toJSON scale),
                      ("log_alpha", J.toJSON logAlpha), ("alpha_tag", J.toJSON (tag alpha)),
                      ("alpha", number alpha), ("latent_log_density", number basePrior),
                      ("latent_log_density_tag", J.toJSON (tag basePrior)),
                      ("derived_log_density", number derivedDensity),
                      ("derived_log_density_tag", J.toJSON (tag derivedDensity)),
                      ("gamma_rates", J.toJSON rates)]
  T.putStrLn (J.cjsonToText (J.toCJSON row))

grid :: [Double]
grid = [-20, -10, -6, -3, 0, 6, 20, 50, 100, 500, 700, 709, 709.7, 709.78, 709.79, 710, 750, 1000, 5000]

main = mapM_ (\prior -> mapM_ (probe prior) grid)
  [("package",6,2), ("centered",0,1), ("broad",0,2)]
