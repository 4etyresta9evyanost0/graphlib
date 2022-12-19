using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Windows;
using System.Threading;
using System.Collections;

using System.Windows.Controls;

namespace Graphs
{
    public class Vertex
    {
        public string Name { get; set; }
        public int Index { get; set; }

        //public Point Location { get; set; }

        public Vertex(string name, int index)
        {
            Name = name;
            Index = index;
        }

        public double Weight { get; set; }

        public Graph Graph { get; protected internal set; }

        public System.Windows.Controls.Grid Grid { get; set; }
        public System.Windows.Shapes.Ellipse Ellipse { get; set; }

        public Edge[] Begins
        {
            get
            {
                return Graph.Edges.FindAll((x) => x.Begin == this).ToArray();
            }
        }

        public Edge[] Ends
        {
            get
            {
                return Graph.Edges.FindAll((x) => x.End == this).ToArray();
            }
        }


        public Edge[] Edges
        {
            get
            {
                var arr = new HashSet<Edge>(Begins.Concat(Ends));
                //var dist = arr.Distinct().ToArray();
                //arr.ExceptWith(dist);
                return arr.ToArray();
            }
        }

        //public Edge[] Edges
        //{
        //    get => Graph.Edges.FindAll((x) => x.Begin == this || x.End == this).ToArray();
        //}

        public int[] AdjacencyRow
        {

            get
            {
                int[] row = new int[Graph.Vertices.Count];

                Begins.ToList().ForEach((x) => row[x.End.Index] = 1);
                Ends.ToList().ForEach((x) => {
                    if (!x.HasDirection)
                        row[x.Begin.Index] = 1;
                });

                return row;
                //int[] row = new int[Graph.Vertices.Count];
                //for (int i = 0; i < Graph.Vertices.Count; i++)
                //{
                //    if (Graph.Vertices[i].Begins.ToList(). )
                //    {

                //        row[Graph.Edges[i].End.Index] = 1;
                //    }
                //}
                //return row;
            }
        }
    }

    public class Edge : IComparable
    {
        public string Name { get; set; }
        public int Index { get; set; }
        public bool HasDirection { get; set; }
        public bool HasWeight { get { return Weight != null; } }
        public double? Weight { get; set; }

        public bool isLoop
        {
            get { return Begin == End; }
        }

        public Vertex Begin { get; set; }
        public Vertex End { get; set; }

        public System.Windows.Controls.Grid Grid { get; set; }
        public System.Windows.Shapes.Line Line { get; set; }

        public Vertex GetOpposite(Vertex v)
        {
            if (v == Begin)
                return End;
            else if (v == End)
                return Begin;
            else
                return null;
                
        }

        public Edge(string name, int index, bool hasDirection)
        {
            Name = name;
            Index = index;
            HasDirection = hasDirection;
        }

        public Edge(string name, int index, bool hasDirection, Vertex oneVertex)
        {
            Name = name;
            Index = index;
            HasDirection = hasDirection;
            Begin = oneVertex;
            End = oneVertex;
        }

        public Edge(string name, int index, bool hasDirection, Vertex begin, Vertex end)
        {
            Name = name;
            Index = index;
            HasDirection = hasDirection;
            Begin = begin;
            End = end;
        }

        public Edge(string name, int index, bool hasDirection, Vertex begin, Vertex end, double? weight)
        {
            Name = name;
            Index = index;
            HasDirection = hasDirection;
            Begin = begin;
            End = end;
            Weight = weight;
        }

        public int CompareTo(object value)
        {
            if (value == null)
            {
                return 1;
            }
            //if (value is Edge num)
            if ((Edge)value != null )
            {
                return CompareTo((Edge)value);
            }
            throw new ArgumentException("Must be edge");
        }

        public int CompareTo(Edge value)
        {
            if (value == null)
            {
                return 1;
            }
            if (Weight < value.Weight || (!HasWeight && value.HasWeight))
            {
                return -1;
            }
            if (Weight > value.Weight || (HasWeight && !value.HasWeight))
            {
                return 1;
            }
            return 0;
        }

        public Graph Graph { get; protected internal set; }
    }

    struct SimpleEdge
    {
        public string Начало { get; set; }
        public string Конец { get; set; }
        public double? Вес { get; set; }

        public SimpleEdge(string b, string e, double? w)
        {
            Начало = b;
            Конец = e;
            Вес = w;
        }
    }


    public enum GraphMatrix { Adjacency = 0, Incidence = 1 }
    public enum GraphType { Empty = 0, Undirected = 1, Directed = 2, Mixed = 3 }

    public class Graph
    {
        public double AllWeights
        {
            get
            {
                double sum = 0;
                if (Edges.Any((x) => !x.HasWeight))
                    return -1;
                Edges.ForEach((x) => sum += x.Weight ?? 0);
                return sum;
            }
        }
        public List<Vertex> Vertices { get; set; } = new List<Vertex>();
        public List<Edge> Edges { get; set; } = new List<Edge>();
        public GraphType Type
        {
            get
            {
                if (Edges.Count == 0 || Vertices.Count == 0)
                    return GraphType.Empty;

                GraphType gt = (Edges[0].HasDirection ? GraphType.Directed : GraphType.Undirected);
                for (int i = 1; i < Edges.Count; i++)
                {
                    if (gt != (Edges[i].HasDirection ? GraphType.Directed : GraphType.Undirected))
                        return GraphType.Mixed;
                }
                return gt;
            }
        }

        public Graph() { }

        public Graph(Edge[] edges, Vertex[] vertices)
        {
            Edges = edges.ToList();
            Vertices = vertices.ToList();

            for (int i = 0; i < Vertices.Count; i++)
            {
                Vertices[i].Graph = this;
                Vertices[i].Index = i;
            }

            for (int i = 0; i < Edges.Count; i++)
            {
                Edges[i].Graph = this;
                Edges[i].Index = i;
            }
        }

        public Graph(int?[,] mat, GraphType type) : base()
        {
            for (int i = 0; i < mat.GetLength(0); i++)
            {
                Vertex e = new Vertex($"X{i + 1}", i);
                Vertices.Add(e);
            }
            if (mat.GetLength(0) == mat.GetLength(1))
            {   // По смежной матрице

                int c = 0;
                switch (type)
                {
                    case GraphType.Directed:
                        for (int i = 0; i < mat.GetLength(0); i++)
                        {
                            for (int j = 0; j < mat.GetLength(1); j++)
                            {
                                if (mat[i, j] != 0)
                                    Edges.Add(new Edge($"U{++c}", c, true, Vertices[i], Vertices[j], mat[i, j]));
                            }
                        }
                        break;
                    case GraphType.Undirected:
                        for (int i = 0; i < mat.GetLength(0); i++)
                        {
                            for (int j = i; j < mat.GetLength(1); j++)
                            {
                                if (mat[i, j] != 0)
                                    Edges.Add(new Edge($"U{++c}", c, false, Vertices[i], Vertices[j], mat[i, j]));
                            }
                        }
                        break;

                    default:
                        throw new NotImplementedException("Нельзя использовать другие типы матриц.");
                }
            }
            else
            {   // По матрице инцидентности
                for (int i = 0; i < mat.GetLength(1); i++)
                {
                    int[] ints = new int[mat.GetLength(0)];
                    for (int j = 0; j < ints.Length; j++)
                    {
                        ints[j] = mat[j, i] ?? 0;
                    }

                    var begin = ints.ToList().FindIndex((x) => x == 1);
                    var end = ints.ToList().FindLastIndex((x) => (type == GraphType.Undirected ? x == 1 : x == -1));

                    if (end != -1)
                        Edges.Add(new Edge($"U{i + 1}", i, Convert.ToBoolean((int)type - 1), Vertices[begin], Vertices[end]));
                    else
                        Edges.Add(new Edge($"U{i + 1}", i, Convert.ToBoolean((int)type - 1), Vertices[begin]));
                }
            }

            // Инициализация
            Vertices.ForEach((x) => x.Graph = this);
            Edges.ForEach((x) => x.Graph = this);
        }

        public int[,] AdjacencyMatrix
        {
            get
            {
                int[,] mat = new int[Vertices.Count, Vertices.Count];

                for (int i = 0; i < Vertices.Count; i++)
                {
                    for (int j = 0; j < Vertices.Count; j++)
                    {
                        mat[i, j] = Vertices[i].AdjacencyRow[j];
                    }
                }

                return mat;
            }
        }

        public int[,] IncidenceMatrix
        {
            get
            {
                int[,] mat = new int[Vertices.Count, Edges.Count];

                for (int i = 0; i < Edges.Count; i++)
                {
                    mat[Edges[i].Begin.Index, i] = 1;
                    mat[Edges[i].End.Index, i] = (Edges[i].HasDirection && Edges[i].Begin.Index != Edges[i].End.Index ? -1 : 1);
                }

                return mat;
            }
        }

        public int[,] ReachabilityMatrix
        {
            get
            {
                if (Type != GraphType.Directed)
                {
                    return null;
                }
                int[,] mat = new int[Vertices.Count, Vertices.Count];
                for (int i = 0; i < Vertices.Count; i++)
                    for (int j = 0; j < Vertices.Count; j++)
                        mat[i, j] = AdjacencyMatrix[i, j];


                for (int k = 0; k < Vertices.Count; k++)
                    for (int i = 0; i < Vertices.Count; i++)
                        for (int j = 0; j < Vertices.Count; j++)
                            mat[i, j] = mat[i, j] | (mat[i, k] & mat[k, j]);

                for (int i = 0; i < Vertices.Count; i++)
                    mat[i, i] = 1;

                return mat;
            }
        }

        public string AdjacencyMatrixStr
        {
            get { return CommonFunctions.GetStrMatrix(AdjacencyMatrix); }
        }

        public string IncidenceMatrixStr
        {
            get { return CommonFunctions.GetStrMatrix(IncidenceMatrix); }
        }

        public string ReachabilityMatrixStr
        {

            get { return CommonFunctions.GetStrMatrix(ReachabilityMatrix); }
        }

        enum V { W = 1, G = 2, B = 3 }

        bool DFS(Vertex v, List<Vertex> visitedV, List<Edge> visitedE)
        {
            if (visitedV.Contains(v))
                return true;
            visitedV.Add(v);
            for (int i = 0; i < v.Edges.Length; i++)
            {
                if (visitedE.Contains(v.Edges[i]))
                    continue;
                visitedE.Add(v.Edges[i]);
                var v2 = v.Edges[i].GetOpposite(v);
                if (DFS(v2, visitedV, visitedE))
                    return true;
            }

            return false;
        }

        public bool HasLoops
        {
            get
            {
                List<Vertex> visitedV = new List<Vertex>();
                List<Edge> visitedE = new List<Edge>();


                Vertex ve;
                //bool b = false;
                while ((ve = Vertices.FirstOrDefault((x) => !visitedV.Contains(x))) != null)
                {
                    if (DFS(ve, visitedV, visitedE))
                        return true;
                }
                        
                return false;
            }
        }

        void EndingFunc(MatrixGraphsWPF.CanvasWindow ue, int timeInTheEnd)
        {
            Thread.Sleep(timeInTheEnd);

            foreach (var v in Vertices)
                ue.Dispatcher.Invoke(() => v.Ellipse.Fill = System.Windows.Media.Brushes.SlateGray);

            foreach (var e in Edges)
                ue.Dispatcher.Invoke(() => e.Line.Stroke = System.Windows.Media.Brushes.SteelBlue);

            Thread.Sleep(timeInTheEnd / 10);

            foreach (var v in Vertices)
                ue.Dispatcher.Invoke(() => v.Ellipse.Fill = System.Windows.Media.Brushes.White);

            foreach (var e in Edges)
                ue.Dispatcher.Invoke(() => e.Line.Stroke = System.Windows.Media.Brushes.DodgerBlue);
        }

        bool GrDFS(Vertex v, MatrixGraphsWPF.CanvasWindow ue, int timeBetweenRecursions, int timeInTheEnd,
            List<Vertex> visitedV, List<Edge> visitedE)
        {
            Thread.Sleep(timeBetweenRecursions);
            foreach (var ve in visitedV)
                ue.Dispatcher.Invoke(() => ve.Ellipse.Fill = System.Windows.Media.Brushes.Black);
            ue.Dispatcher.Invoke(() => v.Ellipse.Fill = System.Windows.Media.Brushes.Gray);
            if (visitedV.Contains(v))
            {
                ue.Dispatcher.Invoke(() => v.Ellipse.Fill = System.Windows.Media.Brushes.Yellow);
                return true;
            }

            visitedV.Add(v);
            for (int i = 0; i < v.Edges.Length; i++)
            {
                foreach (var ed in visitedE)
                    ue.Dispatcher.Invoke(() => ed.Line.Stroke = System.Windows.Media.Brushes.Black);
                var e = v.Edges[i];
                ue.Dispatcher.Invoke(() => e.Line.Stroke = System.Windows.Media.Brushes.Gray);
                if (visitedE.Contains(e))
                {
                    ue.Dispatcher.Invoke(() => e.Line.Stroke = System.Windows.Media.Brushes.Yellow);
                    continue;
                }

                visitedE.Add(e);
                var v2 = e.GetOpposite(v);
                if (GrDFS(v2, ue, timeBetweenRecursions, timeInTheEnd, visitedV, visitedE))
                {
                    return true;
                }
            }

            return false;
        }

        public bool GraphicalDFS(MatrixGraphsWPF.CanvasWindow ue, int timeBetweenRecursions, int timeInTheEnd, UIElement[] toInvisible = null)
        {
            foreach (var v in Vertices)
                ue.Dispatcher.Invoke(() => v.Ellipse.Fill = System.Windows.Media.Brushes.White);

            if (toInvisible != null)
                foreach (var el in toInvisible)
                    ue.Dispatcher.Invoke(() => el.Visibility = Visibility.Hidden);

            //if (toInvisible != null)
            foreach (var e in Edges)
                ue.Dispatcher.Invoke(() => e.Line.Visibility = e.Grid.Visibility = Visibility.Visible);


            foreach (var e in Edges)
                ue.Dispatcher.Invoke(() => e.Line.Stroke = System.Windows.Media.Brushes.DodgerBlue);



            List<Vertex> visitedV = new List<Vertex>();
            List<Edge> visitedE = new List<Edge>();    

            Vertex ver;
            while ((ver = Vertices.FirstOrDefault((x) => !visitedV.Contains(x))) != null)
            {
                if (GrDFS(ver, ue, timeBetweenRecursions, timeInTheEnd, visitedV, visitedE))
                {
                    EndingFunc(ue, timeInTheEnd);
                    return true;
                }
            }

            EndingFunc(ue,timeInTheEnd);
            return false;
        }

        public Graph GetSpanningTreePrim(Vertex begin = null)
        {
            if (begin == null)
            {
                begin = Vertices[0];
            }

            if (Type != GraphType.Undirected)
            {
                return null;
            }

            //var ToExecVertices = new List<Vertex>(Vertices.Except(new Vertex[] { begin }));
            var ExecutedVertices = new List<Vertex>();
            var CurrentExecution = new List<Vertex>();
            var nEdges = new List<Edge>();
            

            while(nEdges.Count != Vertices.Count - 1)
            {
                CurrentExecution.Add(begin);
                ExecutedVertices.Add(begin);
                var query = begin.Edges.Where(x => !CurrentExecution.Contains(x.GetOpposite(begin)) && !nEdges.Contains(x));

                if (nEdges.Count == 6) 
                    ;

                if (query.Count() == 0)
                {
                    CurrentExecution.Clear();
                    begin = Vertices.Except(ExecutedVertices).First();
                    continue;
                }

                //var query = from e in begin.Edges
                //        from v in ExecutedVertices
                //        where v != e.GetOpposite(begin)
                //        select e;

                var min = query.Min();

                begin = min.GetOpposite(begin);
                nEdges.Add(min);

                if (ExecutedVertices.Contains(begin))
                {
                    var quer = Vertices.Except(ExecutedVertices);
                    if (quer.Count() == 0)
                        continue;
                    begin = quer.First();
                }
            }

            return new Graph(nEdges.ToArray(),Vertices.ToArray());
        }

        public Graph GetSpanningTreeKruskal()
        {
            if (Type != GraphType.Undirected)
            {
                return null;
            }

            var edges = Edges.ToArray().ToList();
            var vertices = Vertices.ToArray().ToList();

            edges.Sort((x, y) => (int)(x.Weight - y.Weight));

            var nEdges = new List<Edge>();


            for (int i = 0; nEdges.Count < Vertices.Count - 1; i++)
            {
                nEdges.Add(edges[i]);
                if (new Graph(nEdges.ToArray(), Vertices.ToArray()).HasLoops)
                    nEdges.Remove(edges[i]);
            }

            return new Graph(nEdges.ToArray(),Vertices.ToArray());
        }

        public Graph GraphicalGetSpanningTreeKruskal(MatrixGraphsWPF.CanvasWindow ue, int timeBetweenRecursions, int timeInTheEnd)
        {
            if (Type != GraphType.Undirected)
            {
                return null;
            }

            var edges = Edges.ToArray().ToList();

            var vertices = Vertices.ToArray().ToList();

            edges.Sort((x, y) => (int)(x.Weight - y.Weight));

            var nEdges = new List<Edge>();

            for (int i = 0; nEdges.Count < Vertices.Count - 1; i++)
            {
                nEdges.Add(edges[i]);
                var disEdges = edges.Where((x) => !nEdges.Contains(x)).ToList();

                List<UIElement> disableUI = new List<UIElement>();
                foreach (var e in disEdges)
                    disableUI.AddRange(new UIElement[] { e.Line, e.Grid });

                if (new Graph(nEdges.ToArray(), Vertices.ToArray()).GraphicalDFS(ue, timeBetweenRecursions, timeInTheEnd, disableUI.ToArray()))
                    nEdges.Remove(edges[i]);
            }

            return new Graph(nEdges.ToArray(), Vertices.ToArray());
        }


        public Graph GraphicalGetSpanningTreePrim(MatrixGraphsWPF.CanvasWindow ue, TextBox tb, int timeBetweenRecursions,  Vertex begin = null)
        {
            if (begin == null)
            {
                begin = Vertices[0];
            }

            if (Type != GraphType.Undirected)
            {
                return null;
            }

            //var ToExecVertices = new List<Vertex>(Vertices.Except(new Vertex[] { begin }));
            var ExecutedVertices = new List<Vertex>();
            var CurrentExecution = new List<Vertex>();
            var nEdges = new List<Edge>();

            ue.Dispatcher.Invoke(() => tb.Clear());
            while (nEdges.Count != Vertices.Count - 1)
            {
                CurrentExecution.Add(begin);
                ExecutedVertices.Add(begin);
                var query = begin.Edges.Where(x => !CurrentExecution.Contains(x.GetOpposite(begin)) && !nEdges.Contains(x))
;
                ue.Dispatcher.Invoke(() => {
                    tb.Text += begin.Name + ":\r\n";
                    List<Edge> edges = new List<Edge>();
                    

                    foreach (Vertex v in Vertices.Except(CurrentExecution))
                    {
                        string str = "";
                        Edge edge = null;

                        if (v.Edges.Any(x => x.Begin == begin || x.End == begin))
                        {
                            str = "[" + begin.Name + ";" + (edge = v.Edges.First(x => x.Begin == begin || x.End == begin)).Weight + "]";
                            edges.Add(v.Edges.First(x => x.Begin == begin || x.End == begin)); 
                        }
                        else
                        {
                            str = "[0;inf]";
                        }

                        tb.Text += $"  {v.Name} - {str} {(edge != null && query.Min() == edge ? " - min" : "")}\r\n";
                    }
                    tb.Text += "\r\n";
                    }
                );

                if (query.Count() == 0)
                {
                    CurrentExecution.Clear();
                    begin = Vertices.Except(ExecutedVertices).First();
                    continue;
                }

                var min = query.Min();

                begin = min.GetOpposite(begin);
                nEdges.Add(min);
                ue.Dispatcher.Invoke(() => min.Line.Stroke = System.Windows.Media.Brushes.DarkGreen);

                if (ExecutedVertices.Contains(begin))
                {
                    var quer = Vertices.Except(ExecutedVertices);
                    if (quer.Count() == 0)
                        continue;
                    begin = quer.First();
                }

                Thread.Sleep(timeBetweenRecursions);
            }

            return new Graph(nEdges.ToArray(), Vertices.ToArray());
        }

        public double[,] Dijkstra()
        {
            double[,] dij = new double[Vertices.Count, Vertices.Count];
            for (int i = 0; i < Vertices.Count; i++)
            {
                double[] b = DijkstraOneVertex(Vertices[i]);
                for (int j = 0; j < b.Length; j++)
                    dij[i, j] = b[j];
            }

            return dij;
        }

        public double[] DijkstraOneVertex(Vertex s)
        {
            double[] D = new double[Vertices.Count];
            bool[] used = new bool[Vertices.Count];
            foreach (var v in Vertices)
            {
                D[v.Index] = double.PositiveInfinity;
                used[v.Index] = false;
            }
            D[s.Index] = 0;
            for (int i = 0; i < Vertices.Count; i++)
            {
                Vertex v = null;
                for (int j = 0; j < Vertices.Count; j++)
                {
                    if (!used[j] && (v == null || D[j] < D[v.Index]))
                        v = Vertices[j];
                }
                if (double.IsInfinity(D[v.Index]))
                    break;
                used[v.Index] = true;
                foreach (var e in v.Edges) // релаксация по всем рёбрам из v;
                {
                    if (D[v.Index] + e.Weight < D[e.GetOpposite(v).Index])
                        D[e.GetOpposite(v).Index] = D[v.Index] + e.Weight ?? 0;
                }
            }

            return D;
        }

        public Tuple<double,double>[] NetworkPlanning()
        {
            //for (int i = 0; i < V.length) 
            // 	if (V[V.length-1][i] != 0) 
            // 		throw new Exception();

            double[] VJobEarly = new double[Vertices.Count];//new double[V.length];
            VJobEarly[0] = 0;
            for (int i = 1; i < Vertices.Count; i++)
            {
                double[] curVJobs = new double[Vertices[i].Ends.Length];
                for (int j = 0; j < curVJobs.Length; j++)
                    curVJobs[j] = VJobEarly[Vertices[i].Ends[j].GetOpposite(Vertices[i]).Index] + (Vertices[i].Ends[j].Weight ?? 0);
                VJobEarly[i] = curVJobs.Max();
            }

            double[] VJobLate = new double[Vertices.Count];
            VJobLate[Vertices.Count - 1] = VJobEarly[Vertices.Count - 1];
            for (int i = Vertices.Count - 2; i >= 0; i--)
            {
                double[] curVJobs = new double[Vertices[i].Begins.Length];
                for (int j = 0; j < curVJobs.Length; j++)
                    curVJobs[j] = VJobLate[Vertices[i].Begins[j].GetOpposite(Vertices[i]).Index] - (Vertices[i].Begins[j].Weight ?? 0);
                VJobLate[i] = curVJobs.Min();
            }

            Tuple<double, double>[] result = new Tuple<double, double>[Vertices.Count];
            for (int i = 0; i < result.Length; i++)
                result[i] = new Tuple<double, double>(VJobEarly[i], VJobLate[i]);

            return result;
        }

        public Tuple<double, double>[] NetworkPlanningGraphical(out string content, out Edge[] arrWay)
        {
            double[] VJobEarly = new double[Vertices.Count];
            VJobEarly[0] = 0;
            content = "Ранние сроки завершения работы:\n\r\n\rtр(1) = 0\n\r";
            for (int i = 1; i < Vertices.Count; i++)
            {
                double[] curVJobs = new double[Vertices[i].Ends.Length];
                content += $"tр({i + 1}) = max{{\n\r";
                for (int j = 0; j < curVJobs.Length; j++)
                {
                    int pastIndex = Vertices[i].Ends[j].GetOpposite(Vertices[i]).Index;
                    double pastValue = VJobEarly[pastIndex];
                    double weightBetween = (Vertices[i].Ends[j].Weight ?? 0);
                    curVJobs[j] = pastValue + weightBetween;
                    content += $"\ttр({pastIndex+1}) + t({pastIndex+1}, {i+1}) = {curVJobs[j]}\n\r";
                }
                VJobEarly[i] = curVJobs.Max();
                content += $"}} = {VJobEarly[i]}\n\r\n\r";
            }

            double[] VJobLate = new double[Vertices.Count];
            VJobLate[Vertices.Count - 1] = VJobEarly[Vertices.Count - 1];

            content += $"\n\rРанние сроки завершения работы:\n\r\n\rtр({Vertices.Count}) = 0\n\r";
            for (int i = Vertices.Count - 2; i >= 0; i--)
            {
                double[] curVJobs = new double[Vertices[i].Begins.Length];
                content += $"tр({i + 1}) = min{{\n\r";
                for (int j = 0; j < curVJobs.Length; j++)
                {
                    int pastIndex = Vertices[i].Begins[j].GetOpposite(Vertices[i]).Index;
                    double pastValue = VJobLate[pastIndex];
                    double weightBetween = (Vertices[i].Begins[j].Weight ?? 0);
                    curVJobs[j] = pastValue - weightBetween;
                    content += $"\ttр({pastIndex + 1}) + t({i + 1},{pastIndex + 1}) = {curVJobs[j]}\n\r";
                }
                VJobLate[i] = curVJobs.Min();
                content += $"}} = {VJobLate[i]}\n\r\n\r";
            }

            Tuple<double, double>[] result = new Tuple<double, double>[Vertices.Count];
            for (int i = 0; i < result.Length; i++)
                result[i] = new Tuple<double, double>(VJobEarly[i], VJobLate[i]);

            var reserve = new double[Vertices.Count];

            content += "\n\rРезерв времени:\n\r\n\r";
            for (int i = 0; i < Vertices.Count; i++)
            {
                content += $"R({i + 1}) = {VJobLate[i]} - {VJobEarly[i]} = {reserve[i] = VJobLate[i] - VJobEarly[i]}\r\n\r\n";
            }

            List<int> minWay = new List<int>();

            for (int i = 0; i < Vertices.Count; i++)
            {
                if (reserve[i] == 0)
                    minWay.Add(i);
            }
            content += "\r\nКритический путь: " + string.Join(" - ", minWay.ConvertAll(x => (++x).ToString()));
            Edge[] edgesMinWay = new Edge[minWay.Count - 1];
            for (int i = 0; i < edgesMinWay.Length; i++)
            {
                edgesMinWay[i] = Vertices[minWay[i]].Begins.Where(n => n.End.Index == minWay[i + 1]).First();
            }

            arrWay = edgesMinWay;
            return result;
        }

        //public double[] DijkstraAlgorithm(Vertex begin = null)
        //{
        //    if (begin == null)
        //        begin = Vertices[0];

        //    List<Edge> execEdges = new List<Edge>();

        //    for (int i = 0; i < Vertices.Count; i++)
        //    {
        //        if (begin == Vertices[i])
        //            continue;

        //        Vertex

        //        execEdges.AddRange(begin.Edges);
        //    }

        //    //Dictionary<Vertex, List<double>> executedVerticies = new Dictionary<Vertex, List<double>>();

        //    //while (executedVerticies.Count != Vertices.Count)
        //    //{
        //    //    List<Tuple<double, Vertex>> neighVerts = new List<Tuple<double, Vertex>>();
        //    //    foreach (var e in begin.Edges)
        //    //    {
        //    //        neighVerts.Add(new Tuple<double,Vertex>(e.Weight, e.GetOpposite(begin)));
        //    //    }



        //    //    //double max = double.MinValue;
        //    //    //int ind = 0;
        //    //    //for (int i = 0; i < verts.Count; i++)
        //    //    //    if (max < verts[i].Item1)
        //    //    //    {
        //    //    //        ind = i;
        //    //    //        max = verts[i].Item1;
        //    //    //    }


        //    //}



        //    return null;
        //}

    }



}
