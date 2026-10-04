{-# LANGUAGE OverloadedStrings #-}
module Main where
import qualified Data.JSON as J
import qualified Data.Text as Text
import qualified Data.Text.IO as T
import SModel.ASRV (gammaRatesMean)
import Probability.Distribution.Discrete (values)

probe line = do
  let fields = words line
      ident = fields !! 0
      alpha = if fields !! 1 == "Infinity" then 1/0 else read (fields !! 1) :: Double
      raw = values (gammaRatesMean alpha 4)
      adjusted = map (/ (sum raw / 4)) raw
      checks = map (\x -> x == (read (Text.unpack (J.cjsonToText (J.toCJSON x))) :: Double)) adjusted
      result = J.Object [("id", J.toJSON ident), ("alpha_token", J.toJSON (fields !! 1)), ("alpha", J.toJSON alpha),
                         ("raw_rates", J.toJSON raw), ("model_rates", J.toJSON adjusted),
                         ("roundtrips", J.toJSON checks),
                         ("original_encoding", J.toJSON (Text.unpack (J.encode adjusted)))]
  T.putStrLn (J.cjsonToText (J.toCJSON result))

main = do
  contents <- T.readFile "__INPUT_PATH__"
  mapM_ probe (lines (Text.unpack contents))
