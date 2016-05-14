/*
Solving this problem.

https://www.hackerrank.com/challenges/castle-on-the-grid

A graph is composed of vertices connected by edges. Circles connected by lines.

A road map can be represented by graph.

A maze can be represented by graph.

A course pre-reqs can be represented by graph.

A tree is a special case of a graph where vertices only have 1 predecessor.

To implement a graph we need

A struct for edges. edges have an end vertex and a cost.
A struct for vertices. Vertices have a label and a list of edges. + isvisited, previousVertex.

We need structs for adding vertices to the graph. And connecting vertices with edges.

To implement get cheapest path we need a priority queue.

*/


#include <stdio.h>
#include <stdlib.h>

// we need a list to use as queue in the shortest path algo.

struct _LIST_NODE {
  struct _LIST_NODE *Next;
  int               Priority;
  void              *Data;
};

typedef struct _LIST_NODE LIST_NODE;

struct _LIST {
  LIST_NODE *Head;
  LIST_NODE *Tail;
};

typedef struct _LIST LIST;

// creates a list node.
LIST_NODE *
NewListNode(void *Data, int Priority) {
  LIST_NODE *Node;

  Node = calloc(1, sizeof(LIST_NODE));

  if(Node == NULL) {
    printf("NewListNode. out of mem.\n");
    return NULL;
  }

  Node->Priority = Priority;
  Node->Data = Data;

  return Node;
}


// Appends a node to the tail of the list.
int
ListTailAdd(LIST *List, LIST_NODE *Node) {

  if(List == NULL || Node == NULL) {
    printf("ListAppend. Invalid arg.\n");
    return 1;
  }

  Node->Next = NULL;

  if(List->Head == NULL) { // empty list add.
    List->Head = Node;
    List->Tail = Node;
    return 0;
  }

  // tail add.
  List->Tail->Next = Node;
  List->Tail = Node;

  return 0;
}



// remove head of list
LIST_NODE *
ListHeadRemove(LIST *List) {
  LIST_NODE *Node;

  if(List == NULL) {
    printf("ListHeadRemove. Invalid arg.\n");
    return NULL;
  }

  if(List->Head == NULL) {
    return NULL;
  }

  Node = List->Head;

  List->Head = Node->Next;
  Node->Next = NULL;

  return Node;
}

// remove tail of list
LIST_NODE *
ListTailRemove(LIST *List) {
  LIST_NODE *Node, *OldTail;

  if(List == NULL) {
    printf("ListTailRemove. Invalid arg.\n");
    return NULL;
  }

  if(List->Head == NULL) { // empty list
    return NULL;
  }

  OldTail = List->Tail;

  if(List->Head == List->Tail) { // 1 node in list
    List->Tail = NULL;
    List->Head = NULL;
    return OldTail;
  }

  // Note: to use the list as a stack with O(1) pop() we need a doubly linked list.

  // Find the node before the tail, this will be the new tail
  for(Node = List->Head; Node != NULL && Node->Next != List->Tail; Node = Node->Next)
    ;

  Node->Next = NULL;
  List->Tail = Node;

  return OldTail;
}

int
IsEmpty(LIST *List){
  if(List == NULL) {
    printf("IsEmpty. Invalid arg.\n");
    return 1;
  }

  return List->Head == NULL;
}

// Inserts a node in a sorted list at position determined by Node->Priority.
// The list must be sorted from lowest to highest.
int
ListInsertSorted(LIST *List, LIST_NODE *Node) {
  LIST_NODE *NodeBefore;
  LIST_NODE *NodeAfter;

  if(List == NULL || Node == NULL) {
    printf("ListAppend. Invalid arg.\n");
    return 1;
  }

  Node->Next = NULL;

  if(List->Head == NULL) { // empty list add.
    List->Head = Node;
    List->Tail = Node;
    return 0;
  }

  for(NodeBefore = NULL, NodeAfter = List->Head; NodeAfter != NULL && Node->Priority >= NodeAfter->Priority; NodeBefore = NodeAfter, NodeAfter = NodeAfter->Next){
    ;
  }

  if(NodeBefore == NULL){ // head insert
    Node->Next = List->Head;
    List->Head = Node;
    return 0;
  }

  if(NodeAfter == NULL){ // tail insert
    List->Tail->Next = Node;
    List->Tail = Node;
    return 0;
  }

  // insert between two nodes.
  Node->Next = NodeAfter;
  NodeBefore->Next = Node;

  return 0;
}


// todo: list add sorted. if we need a priority queue to implement cheapest path algo.


// definitions to create a graph.
struct _VERTEX {
  void *Label;
  struct _VERTEX *previousVertex;
  int PathCost;
  int IsVisited;
  struct _VERTEX *Neighbors[8]; // not using edge since unweighted. At most 8 neighbors if moving diagonal is allowed. For now assuming diagonal is not allowed.
  int NumNeighbors;
};

typedef struct _VERTEX VERTEX;

struct _EDGE { // unused.
  VERTEX *EndVertex;
  int    Cost;
};


typedef struct _EDGE EDGE;

struct _GRAPH {
  VERTEX  *Vertices;
  int     NumVertices;
};

typedef struct _GRAPH GRAPH;

struct _POINT { // represents grid points.
  int x;
  int y;
};

typedef struct _POINT POINT;

typedef struct _PRIORITY_QUEUE_ENTRY {
  VERTEX  *PrevVertex;
  int     Cost;
  VERTEX  *Vertex;
} PRIORITY_QUEUE_ENTRY;

PRIORITY_QUEUE_ENTRY *
NewPriorityQueueEntry(VERTEX *PrevVertex, int Cost, VERTEX *Vertex) {
  PRIORITY_QUEUE_ENTRY *Entry;

  Entry = calloc(1, sizeof(PRIORITY_QUEUE_ENTRY));

  if(Entry == NULL) {
    printf("NewPriorityQueueEntry. out of mem.\n");
    return NULL;
  }

  Entry->PrevVertex = PrevVertex;
  Entry->Cost = Cost;
  Entry->Vertex = Vertex;

  return Entry;
}

// enqueues an entry to a priority queue.
int
PriorityEnqueue(LIST *List, void *Data, int Priority){
  LIST_NODE *Node = NewListNode(Data, Priority);

  if(Node == NULL) {
    printf("PriorityEnqueue. out of mem.\n");
    return 0;
  }

  return ListInsertSorted(List, Node);
}

// enqueues an entry
int
Enqueue(LIST *List, void *Data){
  LIST_NODE *Node = NewListNode(Data, 0);

  if(Node == NULL) {
    printf("Enqueue. out of mem.\n");
    return 0;
  }

  return ListTailAdd(List, Node);
}

// dequeues an entry.
void *
Dequeue(LIST *List){
 LIST_NODE *Node = ListHeadRemove(List);
 void *Data = NULL;

 if(Node != NULL){
  Data = Node->Data;
  free(Node);
 }

 // todo: free Node

 return Data;
}


// push a vertex
int
Push(LIST *List, VERTEX *vertex){
  return Enqueue(List, vertex);
}

// pop a vertex.
VERTEX *
Pop(LIST *List){
 LIST_NODE *Node = ListTailRemove(List);

 return (VERTEX *)Node->Data;
}

// creates a Vertex. unused.
VERTEX *
NewVertex(void *Label) {
  VERTEX *Vertex;

  Vertex = calloc(1, sizeof(VERTEX));

  if(Vertex == NULL) {
    printf("NewVertex. out of mem.\n");
    return NULL;
  }

  Vertex->Label = Label;

  return Vertex;
}
void
PrintPoint(POINT *Point) {
  if(Point == NULL)
    return;

  printf("(%d,%d)",Point->x, Point->y);
}

// Prints the entries in List from head to tail. The data in the list nodes must be a pointer to a PRIORITY_QUEUE_ENTRY.
void
PrintPQEntries(LIST *List)
{
  LIST_NODE            *Node;
  PRIORITY_QUEUE_ENTRY *PQEntry;

  if(IsEmpty(List)){
    printf("Empty List.\n");
    return;
  }

  for(Node = List->Head; Node != NULL; Node = Node->Next) {
    PQEntry = Node->Data;
    printf("[");
    PrintPoint(PQEntry->Vertex->Label);
    printf("|c=%d|p=%d|", PQEntry->Cost, Node->Priority); // todo: check the cost stored in vertex?
    if(PQEntry->PrevVertex != NULL) // PrevVertex can be NULL when the first vertex has no predecessor.
      PrintPoint(PQEntry->PrevVertex->Label);
    else
      printf("(NULL)");
    printf("]->");
  }

  printf("\n");
}

// adds a neighbor (i,j) vertex if (i,j) is not a forbidden cell.
int
AddNeighbor(int *grid,  int N, int i, int j, VERTEX *Vertex, VERTEX *Vertices) {
  int NeighborLabel;

  // invalid neighbor label.
  if(i < 0 || i >= N)
    return 1;

  if(j < 0 || j >= N)
    return 1;



  NeighborLabel = i*N + j;

  if(!grid[NeighborLabel])
    return 1;

  Vertex->Neighbors[Vertex->NumNeighbors++] = &Vertices[NeighborLabel];

  return 0;
}

// creates a graph based on grid.
int
CreateGraphFromGrid(int *grid, int N, POINT  *GridPoints, VERTEX *Vertices, GRAPH  *GridGraph) {
  int i, j, vertexlabel;

  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++) {


        // compute index of vertex from its position.

        vertexlabel = i*N + j; // Identify a vertex by is its index in the input grid.
        GridPoints[vertexlabel].x = i; // A vertex label is the point it represents.
        GridPoints[vertexlabel].y = j;
        Vertices[vertexlabel].Label = &GridPoints[vertexlabel];

      if(grid[vertexlabel]) {
        // connect this vertex with its neighbors. // assume you cant move diagonally.
        if(i-1 >= 0){ //  neighbor below
          AddNeighbor(grid, N, i-1, j, &Vertices[vertexlabel], Vertices);
        }
        if(i+1 < N){ // neighbor above
          AddNeighbor(grid, N, i+1, j, &Vertices[vertexlabel], Vertices);
        }
        if(j-1 >= 0){ // neighbor to left
          AddNeighbor(grid, N, i, j-1, &Vertices[vertexlabel], Vertices);
        }
        if(j+1 < N){ // neighbor to right
          AddNeighbor(grid, N, i, j+1, &Vertices[vertexlabel], Vertices);
        }
      }
    }
  }

  GridGraph->Vertices = Vertices;
  GridGraph->NumVertices = N*N;

  return 0;
}

// for testing
int
PrintNeighbors(VERTEX *Vertex){
  int i;
  POINT *TmpPoint;

  TmpPoint = Vertex->Label;

  printf("(%d,%d)\n",TmpPoint->x, TmpPoint->y);

  for(i = 0; i < Vertex->NumNeighbors; i++){
    TmpPoint = Vertex->Neighbors[i]->Label;

    //printf("(%d,%d)\n", TmpPoint->x, TmpPoint->y);
  }

  return 0;
}

// for testing
int
TestQueue(VERTEX *Vertices, int count) {
  LIST List = {0};
  int i;
  VERTEX *vertex;

  for(i=0; i<count; i++){
    Enqueue(&List, &Vertices[i]);
  }

  while(!IsEmpty(&List)){
    vertex = Dequeue(&List);
    PrintNeighbors(vertex);
  }

  return 0;
}

// test for priority queue.
int
TestPriorityQueue(void){
  LIST List = {0};
  LIST_NODE *Node;
  int i;
  #define ARRAY_SIZE 8
  int a[ARRAY_SIZE] = {0,7,5,34,1,2,1,0};
  int *Data;

  for(i=0; i< ARRAY_SIZE; i++){
    ListInsertSorted(&List, NewListNode(&a[i], a[i]));
  }

  while(!IsEmpty(&List)){
    Data = (int *)Dequeue(&List);
    printf("%d\n", *Data);
  }

  return 0;
}

// compute slope from p0 to p1.
void
GetSlope(POINT *p0, POINT *p1, POINT *slope)
{
  slope->x = p1->x - p0->x;
  slope->y = p1->y - p0->y;
}

// returns 1 if p0 == p1.
int
IsPointsEqual(POINT *p0, POINT *p1)
{
  return p1->x == p0->x && p1->y == p0->y;
}


int
GetNumPathSteps(LIST *pathStack) {
  POINT prevSlope = {0}, currentSlope = {0};
  int numsteps = 0;
  POINT *prevPoint = NULL, *currentPoint = NULL;
  VERTEX *vertex;

  if(pathStack == NULL || IsEmpty(pathStack)){
    return 0;
  }

  // compute the path length given the points along the path. each straight line is 1 step.


  vertex = Pop(pathStack);
  prevPoint = vertex->Label;

  if(IsEmpty(pathStack)){ // 1 vertex along path means 0 steps.
    return 0;
  }

  vertex = Pop(pathStack);
  currentPoint = vertex->Label;
  numsteps = 1; // we have at least two points, this is the first step

  GetSlope(prevPoint, currentPoint, &prevSlope);
  GetSlope(prevPoint, currentPoint, &currentSlope);

  while(prevPoint != NULL && currentPoint != NULL){
    //PrintNeighbors(vertex);

    // if slope changes between the prev pair of points of the current pair, add a step
    if(!IsPointsEqual(&prevSlope, &currentSlope)){
        numsteps++;
    }


    // get the next pair of points.
    prevPoint = currentPoint;

    if(IsEmpty(pathStack)){ // done
      break;
    }

    vertex = Pop(pathStack);
    currentPoint = vertex->Label;

    // store current slope, compute the next one
    prevSlope = currentSlope;
    GetSlope(prevPoint, currentPoint, &currentSlope);
  }

  return numsteps;
}

// returns 1 if slope(a,b) != slope(b,c), 0 otherwise.
int
GetStepCost(VERTEX *a,VERTEX *b,VERTEX *c) { // a = frontVertex.PrevVertex, b = frontVertex, c = some neighbor of frontVertex.
  POINT slope0, slope1;

  if(a == NULL && b == NULL && c != NULL) {
    // should not occur
    printf("GetStepCost() b c null.\n");
    return 0;
  }

  if(a == NULL && b != NULL && c != NULL) {
    // b and c are diff. +1, if equal +0
    return !IsPointsEqual(b->Label, c->Label);
  }

  if(a != NULL && b != NULL && c != NULL) {


    // PrintPoint(a->Label);
    // PrintPoint(b->Label);
    // PrintPoint(c->Label);
    // printf("\n");

    GetSlope(a->Label, b->Label, &slope0);
    GetSlope(b->Label, c->Label, &slope1);
    return !IsPointsEqual(&slope0, &slope1);
  }

  printf("GetStepCost() invalid arg.\n");
  return 0;
}

void
GetShortestPath(GRAPH *graph, int N, POINT *startPoint, POINT *endPoint, LIST *pathStack){
  LIST vertexQueue = {0};
  int done = 0;
  int neighborIndex;
  VERTEX *startVertex, *endVertex, *frontVertex, *nextNeighbor;


  startVertex = &graph->Vertices[startPoint->x*N+startPoint->y];
  endVertex = &graph->Vertices[endPoint->x*N+endPoint->y];

  // mark origin as visited
  startVertex->IsVisited = 1;
  Enqueue(&vertexQueue, startVertex);

  while(!done && !IsEmpty(&vertexQueue)){
    frontVertex = Dequeue(&vertexQueue);
    neighborIndex = 0;
    while(!done && neighborIndex < frontVertex->NumNeighbors){
      nextNeighbor = frontVertex->Neighbors[neighborIndex];

      if(!nextNeighbor->IsVisited){
        nextNeighbor->IsVisited = 1;
        nextNeighbor->PathCost = 1 + frontVertex->PathCost;
        nextNeighbor->previousVertex = frontVertex;
        Enqueue(&vertexQueue, nextNeighbor);
      }

      if(nextNeighbor == endVertex)
        done = 1;

      neighborIndex++;
    }
  }
  // todo: free queue nodes mem.

  // Add each vertex along the path to the pathStack starting from the end to the beginning.

  while(endVertex != NULL){
    //PrintNeighbors(endVertex);

    Push(pathStack, endVertex);
    endVertex = endVertex->previousVertex;
  }

}

void
GetCheapestPath(GRAPH *graph, int N, POINT *startPoint, POINT *endPoint, LIST *pathStack){
  LIST vertexQueue = {0};
  int done = 0, TmpCost = 0;
  int neighborIndex;
  VERTEX *startVertex, *endVertex, *frontVertex, *nextNeighbor;
  PRIORITY_QUEUE_ENTRY *PQEntry;

  startVertex = &graph->Vertices[startPoint->x*N+startPoint->y];
  endVertex = &graph->Vertices[endPoint->x*N+endPoint->y];

  TmpCost = 0;
  PQEntry = NewPriorityQueueEntry(NULL, TmpCost, startVertex);

  PriorityEnqueue(&vertexQueue, PQEntry, TmpCost);

  // todo: add asserts where needed. e.g. for alloc fails.

  while(!done && !IsEmpty(&vertexQueue)){
    PrintPQEntries(&vertexQueue); // debugging.
    PQEntry = Dequeue(&vertexQueue);
    frontVertex = PQEntry->Vertex;

    if(!frontVertex->IsVisited) {
      frontVertex->IsVisited = 1;
      frontVertex->PathCost = PQEntry->Cost;
      frontVertex->previousVertex = PQEntry->PrevVertex;

      free(PQEntry);

      if(frontVertex == endVertex) {
        done = 1;
      } else {

        neighborIndex = 0;

        while(neighborIndex < frontVertex->NumNeighbors) {

          nextNeighbor = frontVertex->Neighbors[neighborIndex];

          TmpCost = GetStepCost(frontVertex->previousVertex, frontVertex, nextNeighbor);

          if(!nextNeighbor->IsVisited){
            TmpCost += frontVertex->PathCost;
            PQEntry = NewPriorityQueueEntry(frontVertex, TmpCost, nextNeighbor);
            PriorityEnqueue(&vertexQueue, PQEntry, TmpCost);
          }

          neighborIndex++;
        }

    }

    }
  }
  // todo: free queue nodes mem.
  printf("endv cost %d\n", endVertex->PathCost);

  // Add each vertex along the path to the pathStack starting from the end to the beginning.

  while(endVertex != NULL){
    //PrintNeighbors(endVertex);

    Push(pathStack, endVertex);
    endVertex = endVertex->previousVertex;
  }

}

// TRUE if TestPoint is among the list of points in Path.
int
IsPointInPath(POINT *TestPoint, LIST *Path){
  LIST_NODE *Node;
  POINT     *Point;
  VERTEX    *Vertex;

  for(Node = Path->Head; Node != NULL; Node = Node->Next){
    Vertex = Node->Data;
    Point = Vertex->Label;
    if(IsPointsEqual(Point, TestPoint)){
      return 1;
    }
  }

  return 0;
}

// Prints the grid and marks all points along the path.
void
PrintPath(int N, int *grid, LIST *Path){
  POINT     Point;
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      Point.x = i;
      Point.y = j;
      if(IsPointInPath(&Point, Path)){
        printf("#");
      } else {
        if(grid[i*N+j])
          printf(".");
        else
          printf("X"); // X
      }

    }
    printf("\n");
  }
}

int
CastleOnGrid(int N, int *grid, POINT *start, POINT *end) {
  LIST pathStack = {0};
  // Plan: Create a graph from the input grid and apply shortest/cheapest path algorithm
  GRAPH  GridGraph = {0};

  // There are N*N points in the grid.
  POINT  *GridPoints = calloc(N*N, sizeof(POINT));

  // Create a vertex for each point.
  VERTEX *Vertices   = calloc(N*N, sizeof(VERTEX));

  // Create the graph.
  CreateGraphFromGrid(grid, N, GridPoints, Vertices, &GridGraph);

  // find shortest path.
  //GetShortestPath(&GridGraph, N, start, end, &pathStack);
  GetCheapestPath(&GridGraph, N, start, end, &pathStack);
  PrintPath(N, grid, &pathStack);

  printf("%d\n", GetNumPathSteps(&pathStack));

  // todo: free mem.

  return 0;
}

// next : test case 3 fails- might need to use cheapest path.
/**
Problem, shortest path algo. does not work in all cases. There can be 2 distinct paths
with the same length but a different number of steps.

You might think using the cheapest path algo. is the solution but consider that
the definition of cost in this problem is different from a usual weighted graph
where the cost of a move from one vertex to another depends only on the value
stored in the edge that connects them. To make a choice of next vertex to visit
you simply sort neighbors by their edge cost.

In this case, the cost of an edge from b to c can be 1 or 0 depending on
if the vertex preceding b along the path so far. If b was the first vertex,
then the cost is 1. Otherwise, if b is preceded by a, the cost is 1 if
slope(a,b) != slope(b,c), otherwise the cost is 0.

Idea:

Init. cost of start vertex to 0. Starting prevSlope = NULL.
Before visiting neighbors, assign cost. (We must assign a cost while traversing
due to the definition of)

// returns 1 if slope(a,b) != slope(b,c), 0 otherwise.
int
getStepCost(v a, v b, v c); // a = frontVertex.PrevVertex, b = frontVertex, c = some neighbor of frontVertex.

// assign a cost to each neighbor,
void
SetNeighborCost(v front);

void
SortByCost(v neighbors, numneighbors);


// testcase 2 answer = 13, my computation is 15.

###############################

Debugging...

Theories.

-We need a modified cheapest path algo. since the cost of next edge is not fixed, it depends on the previous edge along the current path.

-There is a bug in my cheapest path algo.

-My step count logic is off.

-Found a way to shave off one step in test5.txt. Go down from origin.

-Next: step through cheapest path algo. on paper to spot bug.
We never go thrugh 1,2 2,2 3,2 for test6
**/

int
main(void) {
  char lines[100][102]; // 102 = 100 + \n + \0
  int N = 3;

  int grid[10000];

  POINT start = {0,0}, end={0,2};

  //TestPriorityQueue();

  scanf("%d\n",&N);

  for(int i = 0; i < N; i++){
    fgets(&lines[i][0], 102, stdin);
  }

  scanf("%d %d %d %d", &start.x, &start.y, &end.x, &end.y);

  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      grid[i*N+j] = (lines[i][j] == '.');
    }
  }
  CastleOnGrid(N, grid, &start, &end);

  return 0;
}
