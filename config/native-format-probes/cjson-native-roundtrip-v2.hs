module Main where
import qualified Data.JSON as J
import qualified Data.Text.IO as T
import qualified Data.Text as Text
values :: [Double]
values = [2.34e-10, 2.34e-11, 2.34e-20, 2.34e-50, 2.34e-100, 2.34e-101, 2.34e10, 2.34e20, 0.0, 1.0, 0.25, 4.0]
checks :: [Bool]
checks = map (\x -> x == (read (Text.unpack (J.cjsonToText (J.toCJSON x))) :: Double)) values
main = do { T.putStrLn (J.cjsonToText (J.toCJSON values)); T.putStrLn (J.cjsonToText (J.toCJSON checks)) }
