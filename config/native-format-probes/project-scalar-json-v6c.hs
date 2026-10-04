-- BEGIN project scalar JSON logger v6
-- These functions inspect already logged Values, never sample or change the model.
projectV6QualityKey = Text.pack "__project_scalar_v6_quality__"

projectV6Leaves path (J.Object fields) = concat
  [projectV6Leaves (path ++ [J.String key]) value | (JSONInternal.Key key, value) <- fields]
projectV6Leaves path (J.Array values) = concat
  [projectV6Leaves (path ++ [J.INumber index]) value | (index, value) <- zip [0..] values]
projectV6Leaves path value = [(path, value)]

projectV6NonfiniteKind x
  | isNaN x = Just (Text.pack "nan")
  | isInfinite x && x > 0 = Just (Text.pack "positive_infinity")
  | isInfinite x = Just (Text.pack "negative_infinity")
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

projectV6Special kind = J.String (Text.append (Text.pack "__project_scalar_v6__:") kind)

projectV6Sanitize (J.Object fields) = J.Object [(key, projectV6Sanitize value) | (key, value) <- fields]
projectV6Sanitize (J.Array values) = J.Array (map projectV6Sanitize values)
projectV6Sanitize (J.FNumber x) = case projectV6NonfiniteKind x of
  Just kind -> projectV6Special kind
  Nothing -> J.FNumber x
projectV6Sanitize J.Null = projectV6Special (Text.pack "literal_null")
projectV6Sanitize (J.String text) =
  if take 22 (Text.unpack text) == "__project_scalar_v6__:"
  then error "Reserved project scalar v6 string collision"
  else J.String text
projectV6Sanitize value = value

projectV6ContextValue fields =
  if any (\(JSONInternal.Key key, _) -> key == projectV6QualityKey) fields
  then error "Reserved project scalar v6 context field collision"
  else J.toCJSON (projectV6Sanitize (J.Object (fields ++ [(J.toJSONKey projectV6QualityKey, projectV6Quality fields)])))

projectV6EvaluateContext action _ context = do
  fields <- runContextAction action context
  return (projectV6ContextValue fields)

projectV6EncodeRecord iteration contextFields parameters = J.pairs
  ((J.toJSONKey "iter") .= iteration <>
   E.pair (J.toJSONKey "statistics//") (J.cjsonToEncoding contextFields) <>
   E.pair (J.toJSONKey "parameters//") (J.cjsonToEncoding (J.toCJSON (projectV6Sanitize (J.Object parameters)))) <>
   E.pair (J.toJSONKey "numericParameterQuality//") (J.cjsonToEncoding (J.toCJSON (projectV6Quality parameters))))

projectV6LogSample handle parameters iteration contextFields = do
  T.hPutStrLn handle (J.fromEncoding (projectV6EncodeRecord iteration contextFields parameters))
  hFlush handle

projectScalarJSONLoggerV6 filename = do
  handle <- openFile filename WriteMode
  hPutStrLn handle "{\"fields\":[\"iter\",\"prior\",\"likelihood\",\"posterior\"],\"nested\":true,\"format\":\"MCON\",\"version\":\"0.2\",\"projectScalarSchema\":\"native-cjson-explicit-special-values-v6\"}"
  hFlush handle
  return (\loggerValues ->
    (projectV6EvaluateContext (contextLogValues loggerValues),
     projectV6LogSample handle (parameterLogValues loggerValues)))
-- END project scalar JSON logger v6
