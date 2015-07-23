var flare = 
{
 name: "flare",
 children: [
  {
   name: "analytics",
   children: [
    {
     name: "cluster",
     children: [
      {name: "AgglomerativeCluster", size: 3938},
      {name: "CommunityStructure", size: 3812},
      {name: "HierarchicalCluster", size: 6714},
      {name: "MergeEdge", size: 743}
     ]
    },
    {
     name: "graph",
     children: [
      {name: "BetweennessCentrality", size: 3534},
      {name: "LinkDistance", size: 5731},
      {name: "MaxFlowMinCut", size: 7840},
      {name: "ShortestPaths", size: 5914},
      {name: "SpanningTree", size: 3416}
     ]
    },
    {
     name: "optimization",
     children: [
      {name: "AspectRatioBanker", size: 7074}
     ]
    }
   ]
  },
  {
   name: "animate",
   children: [
    {name: "Easing", size: 17010},
    {name: "FunctionSequence", size: 5842},
    {
     name: "interpolate",
     children: [
      {name: "ArrayInterpolator", size: 1983},
      {name: "ColorInterpolator", size: 2047},
      {name: "DateInterpolator", size: 1375},
      {name: "Interpolator", size: 8746},
      {name: "MatrixInterpolator", size: 2202},
      {name: "NumberInterpolator", size: 1382},
      {name: "ObjectInterpolator", size: 1629},
      {name: "PointInterpolator", size: 1675},
      {name: "RectangleInterpolator", size: 2042}
     ]
    },
    {name: "ISchedulable", size: 1041},
    {name: "Parallel", size: 5176},
    {name: "Pause", size: 449},
    {name: "Scheduler", size: 5593},
    {name: "Sequence", size: 5534},
    {name: "Transition", size: 9201},
    {name: "Transitioner", size: 19975},
    {name: "TransitionEvent", size: 1116},
    {name: "Tween", size: 6006}
   ]
  },
  {
   name: "data",
   children: [
    {
     name: "converters",
     children: [
      {name: "Converters", size: 721},
      {name: "DelimitedTextConverter", size: 4294},
      {name: "GraphMLConverter", size: 9800},
      {name: "IDataConverter", size: 1314},
      {name: "JSONConverter", size: 2220}
     ]
    },
    {name: "DataField", size: 1759},
    {name: "DataSchema", size: 2165},
    {name: "DataSet", size: 586},
    {name: "DataSource", size: 3331},
    {name: "DataTable", size: 772},
    {name: "DataUtil", size: 3322}
   ]
  },
  {
   name: "display",
   children: [
    {name: "DirtySprite", size: 8833},
    {name: "LineSprite", size: 1732},
    {name: "RectSprite", size: 3623},
    {name: "TextSprite", size: 10066}
   ]
  },
  {
   name: "flex",
   children: [
    {name: "FlareVis", size: 4116}
   ]
  },
  {
   name: "physics",
   children: [
    {name: "DragForce", size: 1082},
    {name: "GravityForce", size: 1336},
    {name: "IForce", size: 319},
    {name: "NBodyForce", size: 10498},
    {name: "Particle", size: 2822},
    {name: "Simulation", size: 9983},
    {name: "Spring", size: 2213},
    {name: "SpringForce", size: 1681}
   ]
  },
  {
   name: "query",
   children: [
    {name: "AggregateExpression", size: 1616},
    {name: "And", size: 1027},
    {name: "Arithmetic", size: 3891},
    {name: "Average", size: 891},
    {name: "BinaryExpression", size: 2893},
    {name: "Comparison", size: 5103},
    {name: "CompositeExpression", size: 3677},
    {name: "Count", size: 781},
    {name: "DateUtil", size: 4141},
    {name: "Distinct", size: 933},
    {name: "Expression", size: 5130},
    {name: "ExpressionIterator", size: 3617},
    {name: "Fn", size: 3240},
    {name: "If", size: 2732},
    {name: "IsA", size: 2039},
    {name: "Literal", size: 1214},
    {name: "Match", size: 3748},
    {name: "Maximum", size: 843},
    {
     name: "methods",
     children: [
      {name: "add", size: 593},
      {name: "and", size: 330},
      {name: "average", size: 287},
      {name: "count", size: 277},
      {name: "distinct", size: 292},
      {name: "div", size: 595},
      {name: "eq", size: 594},
      {name: "fn", size: 460},
      {name: "gt", size: 603},
      {name: "gte", size: 625},
      {name: "iff", size: 748},
      {name: "isa", size: 461},
      {name: "lt", size: 597},
      {name: "lte", size: 619},
      {name: "max", size: 283},
      {name: "min", size: 283},
      {name: "mod", size: 591},
      {name: "mul", size: 603},
      {name: "neq", size: 599},
      {name: "not", size: 386},
      {name: "or", size: 323},
      {name: "orderby", size: 307},
      {name: "range", size: 772},
      {name: "select", size: 296},
      {name: "stddev", size: 363},
      {name: "sub", size: 600},
      {name: "sum", size: 280},
      {name: "update", size: 307},
      {name: "variance", size: 335},
      {name: "where", size: 299},
      {name: "xor", size: 354},
      {name: "_", size: 264}
     ]
    },
    {name: "Minimum", size: 843},
    {name: "Not", size: 1554},
    {name: "Or", size: 970},
    {name: "Query", size: 13896},
    {name: "Range", size: 1594},
    {name: "StringUtil", size: 4130},
    {name: "Sum", size: 791},
    {name: "Variable", size: 1124},
    {name: "Variance", size: 1876},
    {name: "Xor", size: 1101}
   ]
  },
  {
   name: "scale",
   children: [
    {name: "IScaleMap", size: 2105},
    {name: "LinearScale", size: 1316},
    {name: "LogScale", size: 3151},
    {name: "OrdinalScale", size: 3770},
    {name: "QuantileScale", size: 2435},
    {name: "QuantitativeScale", size: 4839},
    {name: "RootScale", size: 1756},
    {name: "Scale", size: 4268},
    {name: "ScaleType", size: 1821},
    {name: "TimeScale", size: 5833}
   ]
  },
  {
   name: "util",
   children: [
    {name: "Arrays", size: 8258},
    {name: "Colors", size: 10001},
    {name: "Dates", size: 8217},
    {name: "Displays", size: 12555},
    {name: "Filter", size: 2324},
    {name: "Geometry", size: 10993},
    {
     name: "heap",
     children: [
      {name: "FibonacciHeap", size: 9354},
      {name: "HeapNode", size: 1233}
     ]
    },
    {name: "IEvaluable", size: 335},
    {name: "IPredicate", size: 383},
    {name: "IValueProxy", size: 874},
    {
     name: "math",
     children: [
      {name: "DenseMatrix", size: 3165},
      {name: "IMatrix", size: 2815},
      {name: "SparseMatrix", size: 3366}
     ]
    },
    {name: "Maths", size: 17705},
    {name: "Orientation", size: 1486},
    {
     name: "palette",
     children: [
      {name: "ColorPalette", size: 6367},
      {name: "Palette", size: 1229},
      {name: "ShapePalette", size: 2059},
      {name: "SizePalette", size: 2291}
     ]
    },
    {name: "Property", size: 5559},
    {name: "Shapes", size: 19118},
    {name: "Sort", size: 6887},
    {name: "Stats", size: 6557},
    {name: "Strings", size: 22026}
   ]
  }
 ]
}