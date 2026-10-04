-- BEGIN project scalar JSON logger v6
-- These functions inspect already logged Values, never sample or change the model.
projectV6QualityKey = Text.pack "__project_scalar_v6_quality__"

projectV6Leaves path (J.Object fields) = concat
  [projectV6Leaves (path ++ [J.String key]) value | (JSONInternal.Key key, value) <- fields]
projectV6Leaves path (J.Array values) = concat
  [projectV6Leaves (path ++ [J.INumber index]) value | (index, value) <- zip [0..] values]
projectV6Leaves path value = [(path, value)]

projectV6NonfiniteKind x
  | x /= x = Just (Text.pack "nan")
  | x == (1 / 0) = Just (Text.pack "positive_infinity")
  | x == ((-1) / 0) = Just (Text.pack "negative_infinity")
  | otherwise = Nothing

projectV6Quality fields = J.Object
  [(J.toJSONKey "numericLeafCount", J.INumber (length numbers)),
   (J.toJSONKey "nonfinite", J.Array nonfinite),
   (J.toJSONKey "literalNullPaths", J.Array nulls)]
  where
    leaves = projectV6Leaves [] (J.Object fields)
    numbers = [path | (path, J.INumber _) <- leaves] ++ [path | (path, J.FNumber _) <- leaves]
    nonfinite = [J.Object [(J.toJSONKey "path", J.Array path),
                          (J.toJSONKey "kind", J.String kind)]
                | (path, J.FNumber x) <- leaves, Just kind <- [projectV6NonfiniteKind x]]
    nulls = [J.Array path | (path, J.Null) <- leaves]

projectV6ContextValue fields =
  if any (\(JSONInternal.Key key, _) -> key == projectV6QualityKey) fields
  then error "Reserved project scalar v6 context field collision"
  else J.toCJSON (J.Object (fields ++ [(J.toJSONKey projectV6QualityKey, projectV6Quality fields)]))

projectV6EvaluateContext action _ context = do
  fields <- runContextAction action context
  return (projectV6ContextValue fields)

projectV6EncodeRecord iteration contextFields parameters = J.pairs
  ((J.toJSONKey "iter") .= iteration <>
   E.pair (J.toJSONKey "statistics//") (J.cjsonToEncoding contextFields) <>
   E.pair (J.toJSONKey "parameters//") (J.cjsonToEncoding (J.toCJSON (J.Object parameters))) <>
   E.pair (J.toJSONKey "numericParameterQuality//") (J.cjsonToEncoding (J.toCJSON (projectV6Quality parameters))))

projectV6LogSample handle parameters iteration contextFields = do
  T.hPutStrLn handle (J.fromEncoding (projectV6EncodeRecord iteration contextFields parameters))
  hFlush handle

projectScalarJSONLoggerV6 filename = do
  handle <- openFile filename WriteMode
  hPutStrLn handle "{\"fields\":[\"iter\",\"prior\",\"likelihood\",\"posterior\"],\"nested\":true,\"format\":\"MCON\",\"version\":\"0.2\",\"projectScalarSchema\":\"native-cjson-explicit-nonfinite-v6\"}"
  hFlush handle
  return (\loggerValues ->
    (projectV6EvaluateContext (contextLogValues loggerValues),
     projectV6LogSample handle (parameterLogValues loggerValues)))
-- END project scalar JSON logger v6
