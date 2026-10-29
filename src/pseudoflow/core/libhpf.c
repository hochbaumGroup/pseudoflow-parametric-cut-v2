/*************************************************************************
 * Hochbaum's Pseudo-flow (HPF) Algorithm for Parametric Minimimum Cut   *
 * ***********************************************************************
 Copyright © 2017. The Regents of the University of California (Regents). All Rights Reserved.

 Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities. Created by Quico Spaen and Dorit S. Hochbaum, Department of Industrial Engineering and Operations Research, University of California, Berkeley. This work is adapted from the Pseudoflow implementation by Bala Chandran and Dorit S. Hochbaum available at https://riot.ieor.berkeley.edu/Applications/Pseudoflow/maxflow.html

 IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

 * The HPF algorithm for finding Minimum-cut in a graph is described in: *
 * [1] D.S. Hochbaum, "The Pseudoflow algorithm: A new algorithm for the *
 * maximum flow problem", Operations Research, 58(4):992-1009,2008.      *
 *                                                                       *
 * The algorithm was found to be fast in theory (see the above paper)    *
 * and in practice (see:                                                 *
 * [2] D.S. Hochbaum and B. Chandran, "A Computational Study of the      *
 * Pseudoflow and Push-relabel Algorithms for the Maximum Flow Problem,  *
 * Operations Research, 57(2):358-376, 2009.                             *
 *                                                                       *
 * and                                                                   *
 *                                                                       *
 * [3] B. Fishbain, D.S. Hochbaum, S. Mueller, "Competitive Analysis of  *
 * Minimum-Cut Maximum Flow Algorithms in Vision Problems,               *
 * arXiv:1007.4531v2 [cs.CV]                                             *
 *                                                                       *
 * The algorithm solves a parametric s-t minimum cut problem. The		 *
 * algorithm finds all breakpoints for which the source set of the		 *
 * minimum cut changes as a function of lambda in the range				 *
 * [lower bound, upper bound] by recursively concluding that the interval  *
 * contains 0, 1, or more breakpoints. If the interval contains more than*
 * 1 breakpoint, then the interval is split into two interval, each of   *
 * which contains at least one breakpoint.								 *
 *                                                                       *
 * Parametric cut/flow problems allow for a linear function with input   *
 * lambda on source or sink adjacent arcs. Arcs that are adjacent to	 *
 * source should be non-decreasing in lambda and sink adjacent arcs		 *
 * should be non-increasing in lambda. The algorithm is able to deal with*
 * the reverse configuration (non-increasing on source adjacent arcs and *
 * non-decreasing on sink adjacent arcs) by flipping source and sink and *
 * reversing the direction of the arcs.									 *
 *                                                                       *
 * Usage:																 *
 * 1. Compile hpf.c with a C-compiler (e.g. gcc)						 *
 * 2. To execute within bash environment:								 *
 *	 <name compiled hpf executable> <path input file> <path output file> *
 *                                                                       *
 * INPUT FILE                                                            *
 * **********                                                            *
 * The input file is assumed to be in a modified DIMACS format:	         *
 * c <comment lines>													 *
 * p <# nodes> <# arcs> <lower bound> <upper bound> <round if negative>  *
 * n <source node> s													 *
 * n <sink node> t														 *
 * a <from-node> <to-node> <constant capacity> <lambda multiplier>		 *
 * where the following conditions are satisfied:						 *
 * - Nodes are labeled 0 .. <# nodes> - 1								 *
 * - <lambda multiplier> is non-negative if <from-node> == <source node> *
 *		and <to-node> != <sink-node>									 *
 * - <lambda multiplier> is non-positive if <from-node> != <source node> *
 *		and <to-node> == <sink-node>									 *
 * - <lambda multiplier> is zero if <from-node> != <source node>		 *
 *		and <to-node> != <sink-node>									 *
 * - <lambda multiplier> can take any value if							 *
 *		<from-node> != <source node> and <to-node> != <sink-node>		 *
 * - <round if negative> takes value 1 if the any negative capacity arc  *
 *		should be rounded to 0, and it takes value 0 otherwise			 *
 *                                                                       *
 * OUTPUT FILE                                                           *
 * ***********                                                           *
 * The solver will generate the following output file:					 *
 * t <time (in sec) read data> <time (in sec) initialize> <time			 *
		(in sec) solve>													 *
 * s <# arc scans> <# mergers> <# pushes> <# relabels > <# gap >		 *
 * p <number of lambda intervals = k>									 *
 * l <lambda upperbound interval 1> ... <lambda upperbound interval k>   *
 * n <node-id> <sourceset indicator intval 1 > .. <indicator intval k>   *
 *                                                                       *
 * Set-up                                                                *
 * ******                                                                *
 * Uncompress the MatlabHPF.zip file into the Matlab's working directory *
 * The zip file contains the following files:                            *
 * hpf.c - source code                                                   *
 * hpf.m - Matlab's help file                                            *
 * hpf.mexmaci - The compiled code for Mac OS 10.0.5 (Intel)/ Matlab     *
 *               7.6.0.324 (R2008a).                                     *
 * hpf.mexw32  - The compiled code for Windows 7 / Matlab 7.11.0.584     *
 *               (R2010b).                                               *
 * demo_general - Short Matlab code that generates small network and     *
 *                computes the minimum flow                              *
 * demo_vision - Short Matlab code that loads a Multiview reconstruction *
 *               vision problem (see: [3]) and computes its minimum cut. *
 * gargoyle-smal.mat - The vision problem.                               *
 *                                                                       *
 * When using this code, please cite:                                    *
 * References [1], [2] and [3] above and:                                *
 * Q. Spaen, B. Fishbain and D.S. Hochbaum, "Hochbaum's Pseudo-flow C	 *
 * Implementation", http://riot.ieor.berkeley.edu/riot/Applications/     *
 * Pseudoflow/maxflow.html                                               *
 *************************************************************************/

#define _CRTDBG_MAP_ALLOC
#include "stdio.h"
#include <assert.h>
//#include <sys/time.h>
//#include <sys/resource.h>
#include "stdlib.h"
#include "time.h"
#include <math.h>
//#include <unistd.h>

/*************************************************************************
Definitions
*************************************************************************/
#define  MAX_LEVELS  300
#define VERSION 3.3

typedef unsigned int uint;
typedef long int lint;
typedef long long int llint;
typedef unsigned long long int ullint;

typedef struct Arc
	{
		struct Node *from;
		struct Node *to;
		double flow;
		double capacity;
		double constant;
		double multiplier;
		uint direction;
	} Arc;

typedef struct Node
	{
		struct Node *parent;
		struct Node *childList;
		struct Node *nextScan;
		Arc **outOfTree;
		Arc *arcToParent;
		struct Node *next;
		struct Node *prev;
		double excess;
		uint visited;
		uint numAdjacent;
		uint number;
		int originalIndex;
		uint label;
		uint numOutOfTree;
		uint nextArc;
	} Node;

typedef struct CutProblem
{
	uint numNodesInList;
  uint numSourceSet;
	uint numSinkSet;
	uint numArcs;
	uint solved;
	double lambdaValue;
	Arc *arcList;
	Node *nodeList;
	double cutValue;
	Node *sourceSet;
	Node *sinkSet;
    char *optimalSourceSetIndicator;
} CutProblem;

typedef struct Root
{
	Node *start;
	Node *end;
} Root;

typedef struct Breakpoint
{
	double lambdaValue;
	//uint* sourceSetIndicator;
	struct Breakpoint *next;
} Breakpoint;

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

/*************************************************************************
Global variables
*************************************************************************/
// tolerance for denominator == 0
static double TOL = 1E-7;
static uint numNodes = 0;
static uint numArcs = 0;
static uint numNodesSuper = 0;
static uint numArcsSuper = 0;
static uint numParametricArcs = 0; // Arcs that go from source or into sink
static uint source;
static uint sourceSuper;
static uint sink;
static uint sinkSuper;
static uint highestStrongLabel = 1;

static uint numBreakpoints = 0;
static uint numArcScans = 0;
static uint numPushes = 0;
static uint numMergers = 0;
static uint numRelabels = 0;
static uint numGaps = 0;

static Node *nodesList = NULL;
static Root *strongRoots = NULL;
static uint *labelCount = NULL;
static Arc *arcList = NULL;
static Arc **paramArcList = NULL;
static Node *nodeListSuper = NULL;
static Arc *arcListSuper = NULL;
static uint lowestPositiveExcessNode = 0;

static double *nodeBreakpoints = NULL;

static Breakpoint *lastBreakpoint = NULL;
static Breakpoint *firstBreakpoint = NULL;

static uint useParametricCut = 1;
static uint roundNegativeCapacity = 0;

static double LAMBDA_LOW;
static double LAMBDA_HIGH;

// memory management
static char * all_sink = NULL;
static char * all_source = NULL;
static int * nodeMap= NULL;
static int * sourceAdjacentArcIndices= NULL;
static int * sinkAdjacentArcIndices= NULL;

static Arc *arcListCache[2]={NULL,NULL};
static Node *nodeListCache[2]={NULL,NULL};;
static Node *sourceSetCache[2]={NULL,NULL};;
static Node *sinkSetCache[2]={NULL,NULL};;
static char *pdifferenceCache = NULL;

static time_t rawtime;
static struct tm *info;


double dabs(double value)
{
	if (value >= 0)
		return value;
	else return -value;
}

double math_min(double x, double y)
{
	if (x >= y)
		return y;
	else return x;
}

double math_max(double x, double y)
{
	if (x >= y)
		return x;
	else return y;
}

void libfree(void* p)
{
	free(p);
}

int isFlow(double flow)
/*************************************************************************
isFlow: We set a threshhold. If the flow value is below the threshhold, we
take it as no flow. Otherwise we take it as a flow
*************************************************************************/
{
	if (flow > 0)
		return 1;
	else return 0;
}

int isExcess(double excess)
/*************************************************************************
isExcess: We set a threshhold. If the absolute value of the excess is within
the threshold, then we take it as nothing. Otherwise we will return the sign
of the excess (deficit if negative).
*************************************************************************/
{
	if (excess < 0)
		return -1;
	else if (excess > 0)
		return 1;
	else return 0;
}


static void createOutOfTree (Node *nd)
{
/*************************************************************************
createOutOfTree
*************************************************************************/
	if (nd->numAdjacent)
	{
		if ((nd->outOfTree = (Arc **) malloc (nd->numAdjacent * sizeof (Arc *))) == NULL)
		{
			printf("Out of memory\n");
			exit(0);
		}
	}
}


static void initializeArc (Arc *ac)
{
/*************************************************************************
initializeArc
*************************************************************************/
	ac->from = NULL;
	ac->to = NULL;
	ac->capacity = 0.0;
	ac->flow = 0.0;
	ac->direction = 1;
	ac->constant = 0.0;
	ac->multiplier = 0.0;
}

static void liftAll (Node *rootNode)
{
/*************************************************************************
liftAll
*************************************************************************/
	Node *temp, *current=rootNode;

	current->nextScan = current->childList;

	-- labelCount[current->label];
	current->label = numNodes;

	for ( ; (current); current = current->parent)
	{
		while (current->nextScan)
		{
			temp = current->nextScan;
			current->nextScan = current->nextScan->next;
			current = temp;
			current->nextScan = current->childList;

			-- labelCount[current->label];
			current->label = numNodes;
		}
	}
}
static void addOutOfTreeNode (Node *n, Arc *out)
{
/*************************************************************************
addOutOfTreeNode
*************************************************************************/
	n->outOfTree[n->numOutOfTree] = out;
	++ n->numOutOfTree;
}

static uint sum_array_uint(uint *array, uint num_elements)
{
    uint sum = 0;
    for (uint i = 0; i< num_elements; i++)
    {
        sum += array[i];
    }
    return sum;
}


static uint sum_array_char(char *array, uint num_elements)
{
    uint sum = 0;
		// printf("Number of elements is %d\n", num_elements);
    for (uint i = 0; i< num_elements; i++)
    {
				// printf("array %d is %c\n", i, array[i]); fflush;
        sum += array[i];
    }
    return sum;
}

static void addToStrongBucket (Node *newRoot, Root *rootBucket)
{
/*************************************************************************
addToStrongBucket
*************************************************************************/
	if (rootBucket->start)
	{
		rootBucket->end->next = newRoot;
		rootBucket->end = newRoot;
		newRoot->next = NULL;
	}
	else
	{
		rootBucket->start = newRoot;
		rootBucket->end = newRoot;
		newRoot->next = NULL;
	}
}

static __inline int addRelationship (Node *newParent, Node *child)
{
/*************************************************************************
addRelationship
*************************************************************************/
	child->parent = newParent;
	child->next = newParent->childList;

#ifdef PREV_CHILDLIST
    if (newParent->childList != NULL)
    {
        newParent->childList->prev = child;
    }
#endif

	newParent->childList = child;
	return 0;
}

static __inline void breakRelationship (Node *oldParent, Node *child)
{
/*************************************************************************
breakRelationship
*************************************************************************/
	Node *current;

	child->parent = NULL;

	if (oldParent->childList == child)
	{
		oldParent->childList = child->next;
		child->next = NULL;
		return;
	}
#ifdef PREV_CHILDLIST
    current = child->prev;
    if (child->next != NULL) child->next->prev = current;

#else
	for (current = oldParent->childList; (current->next != child); current = current->next);
#endif
	current->next = child->next;
	child->next = NULL;
}

static void merge (Node *parent, Node *child, Arc *newArc)
{
/*************************************************************************
merge
*************************************************************************/
	Arc *oldArc;
	Node *current = child, *oldParent, *newParent = parent;

	++ numMergers;

	while (current->parent)
	{
		oldArc = current->arcToParent;
		current->arcToParent = newArc;
		oldParent = current->parent;
		breakRelationship (oldParent, current);
		addRelationship (newParent, current);
		newParent = current;
		current = oldParent;
		newArc = oldArc;
		newArc->direction = 1 - newArc->direction;
	}

	current->arcToParent = newArc;
	addRelationship (newParent, current);
}


static __inline void pushUpward (Arc *currentArc, Node *child, Node *parent, const double resCap)
{
/*************************************************************************
pushUpward
*************************************************************************/
	++ numPushes;

	if (isExcess(resCap-child->excess) >= 0)//(/*(int)*/resCap >= child->excess)
	{
		parent->excess += child->excess;
		currentArc->flow += child->excess;
		child->excess = 0;
		return;
	}

	currentArc->direction = 0;
	parent->excess += resCap;
	child->excess -= resCap;
	currentArc->flow = currentArc->capacity;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, &strongRoots[child->label]);
}


static __inline void pushDownward (Arc *currentArc, Node *child, Node *parent, double flow)
{
/*************************************************************************
pushDownward
*************************************************************************/
	++ numPushes;

	if (isExcess(flow - child->excess) >= 0)//(/*(int)*/flow >= child->excess)
	{
		parent->excess += child->excess;
		currentArc->flow -= child->excess;
		child->excess = 0;
		return;
	}

	currentArc->direction = 1;
	child->excess -= flow;
	parent->excess += flow;
	currentArc->flow = 0;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, &strongRoots[child->label]);
}

static void printCutProblem(CutProblem *p){
    printf("numNodes: %u\n " ,p->numNodesInList);
    printf("numSource %u\n" ,p->numSourceSet);
    printf("numSink: %u\n" ,p->numSinkSet);
    printf("numArcs: %u\n" ,p->numArcs);
    printf("solved: %u\n" ,p->solved);
    printf("lambda:%.12lf\n" ,p->lambdaValue);
    int i;
    printf("[from, to](capacity,constant,multiplier)\n");
    for(i=0;i<p->numArcs;++i)
    {
        printf("[%d,%d](%.12lf,%.12lf,%.12lf)\n",p->arcList[i].from->originalIndex,p->arcList[i].to->originalIndex,p->arcList[i].capacity,p->arcList[i].constant,p->arcList[i].multiplier);
    }
    printf("\n");
    //printArcListInfo(arcList);
    //printNodeListInfo(nodeList);
    printf("%.12lf " ,p->cutValue);
    // printNodeListInfo(sourceSet);
    // printSinkListInfo(sinkSet);
    if (p->solved == 1)
    {
				int sourcenodes = 0;
        for(i=0;i<numNodesSuper;++i)
        {
					sourcenodes += p->optimalSourceSetIndicator[i];
            printf("%u ",p->optimalSourceSetIndicator[i]);
        }
				printf("\n");
				printf("Nodes in source set: %d\n", sourcenodes);
    }
    printf("SourceSuper: %d, SinkSuper: %d\n", sourceSuper, sinkSuper);
    printf("\n");
    printf("\n");
}

static void pushExcess (Node *strongRoot)
{
/*************************************************************************
pushExcess
*************************************************************************/
	Node *current, *parent;
	Arc *arcToParent;
	/*int*/double prevEx=1;

	for (current = strongRoot; (isExcess(current->excess) && current->parent); current = parent)
	{
		parent = current->parent;
		prevEx = parent->excess;

		arcToParent = current->arcToParent;

		if (arcToParent->direction)
		{
			pushUpward (arcToParent, current, parent, (arcToParent->capacity - arcToParent->flow));
		}
		else
		{
			pushDownward (arcToParent, current, parent, arcToParent->flow);
		}
	}

	if ((isExcess(current->excess) > 0) && (isExcess(prevEx) <= 0))
	{
		addToStrongBucket (current, &strongRoots[current->label]);
	}
}


static Arc * findWeakNode (Node *strongNode, Node **weakNode)
{
/*************************************************************************
findWeakNode
*************************************************************************/
	uint i, size;
	Arc *out;

	size = strongNode->numOutOfTree;

	for (i=strongNode->nextArc; i<size; ++i)
	{
		++ numArcScans;
		if (strongNode->outOfTree[i]->to->label == (highestStrongLabel-1))
		{
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->to;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		} else if (strongNode->outOfTree[i]->from->label == (highestStrongLabel-1)) {
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->from;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		}
	}

	strongNode->nextArc = strongNode->numOutOfTree;

	return NULL;
}


static void checkChildren (Node *curNode)
{
/*************************************************************************
checkChildren
*************************************************************************/
	for ( ; (curNode->nextScan); curNode->nextScan = curNode->nextScan->next)
	{
		if (curNode->nextScan->label == curNode->label)
		{
			return;
		}

	}

	-- labelCount[curNode->label];
	++	curNode->label;
	++ labelCount[curNode->label];

	++numRelabels;

	curNode->nextArc = 0;
}


static void simpleInitialization (void)
{
/*************************************************************************
simpleInitialization
*************************************************************************/
	uint i, size;
	Arc *tempArc;

	size = nodesList[source].numOutOfTree;
	for (i=0; i<size; ++i) // Saturating source adjacent nodes
	{
		tempArc = nodesList[source].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->to->excess += tempArc->capacity;
	}

	size = nodesList[sink].numOutOfTree;
	for (i=0; i<size; ++i) // Pushing maximum flow on sink adjacent nodes
	{
		tempArc = nodesList[sink].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->from->excess -= tempArc->capacity;
	}

	nodesList[source].excess = 0; // zeroing source excess
	nodesList[sink].excess = 0;	// zeroing sink excess

	for (i=0; i<numNodes; ++i)
	{
		if (isExcess(nodesList[i].excess) > 0)
		{
		    nodesList[i].label = 1;
			++ labelCount[1];

			addToStrongBucket (&nodesList[i], &strongRoots[1]);
		}
	}

	nodesList[source].label = numNodes;	// Set the source label to n
	nodesList[sink].label = 0;			// set the sink label to 0
	labelCount[0] = (numNodes - 2) - labelCount[1];
}


static Node* getHighestStrongRoot (void)
{
/*************************************************************************
getHighestStrongRoot
*************************************************************************/
	uint i;
	Node *strongRoot;

	for (i=highestStrongLabel; i>0; --i)
	{
		if (strongRoots[i].start)
		{
			highestStrongLabel = i;
			if (labelCount[i-1])
			{
				strongRoot = strongRoots[i].start;
				strongRoots[i].start = strongRoot->next;
				strongRoot->next = NULL;
				return strongRoot;
			}

			while (strongRoots[i].start)
			{
				++ numGaps;

				strongRoot = strongRoots[i].start;
				strongRoots[i].start = strongRoot->next;
				liftAll (strongRoot);
			}
		}
	}

	if (!strongRoots[0].start)
	{
		return NULL;
	}

	while (strongRoots[0].start)
	{
		strongRoot = strongRoots[0].start;
		strongRoots[0].start = strongRoot->next;
		strongRoot->label = 1;
		-- labelCount[0];
		++ labelCount[1];

		++ numRelabels;

		addToStrongBucket (strongRoot, &strongRoots[strongRoot->label]);
	}

	highestStrongLabel = 1;

	strongRoot = strongRoots[1].start;
	strongRoots[1].start = strongRoot->next;
	strongRoot->next = NULL;

	return strongRoot;
}



static void initializeRoot (Root *rt)
{
/*************************************************************************
initializeRoot
*************************************************************************/
	rt->start = NULL;
	rt->end = NULL;
}


static void initializeNode (Node *nd, const uint n)
{
/*************************************************************************
initializeNode
*************************************************************************/
	nd->label = 0;
	nd->excess = 0.0;
	nd->parent = NULL;
	nd->childList = NULL;
	nd->nextScan = NULL;
	nd->nextArc = 0;
	nd->numOutOfTree = 0;
	nd->arcToParent = NULL;
	nd->next = NULL;
	nd->prev= NULL;
	nd->visited = 0;
	nd->numAdjacent = 0;
	nd->number = n;
	nd->originalIndex = -10;
	nd->outOfTree = NULL;
}

static void destroyBreakpoint(Breakpoint *currentBreakpoint)
/*************************************************************************
destroyBreakpoint - Removes breakpoint and subsequent ones
*************************************************************************/
{
	if (currentBreakpoint != NULL)
	{
        /* free sourceset indicator */
    	//free(currentBreakpoint->sourceSetIndicator);
    	//currentBreakpoint->sourceSetIndicator = NULL;

    	/* iterate through breakpoint list */
    	destroyBreakpoint(currentBreakpoint->next);
    	currentBreakpoint->next = NULL;

    	/* free breakpoint */
    	free(currentBreakpoint);
    	currentBreakpoint = NULL;
	}

}

static void freeRoot (Root *rt)
{
/*************************************************************************
freeRoot
*************************************************************************/
	rt->start = NULL;
	rt->end = NULL;
}

static void freeCaches(void)
{

    if(all_sink!=NULL) free(all_sink);
    if(all_source!=NULL) free(all_source);
    if(nodeMap!=NULL) free(nodeMap);
    if(sourceAdjacentArcIndices!=NULL) free(sourceAdjacentArcIndices);
    if(sinkAdjacentArcIndices!=NULL) free(sinkAdjacentArcIndices);
    if(pdifferenceCache!=NULL) free(pdifferenceCache);

    for(int i=0; i<2; i++)
    {
        if(arcListCache[i]!=NULL) free(arcListCache[i]);
        if(nodeListCache[i]!=NULL) free(nodeListCache[i]);
        if(sourceSetCache[i]!=NULL) free(sourceSetCache[i]);
        if(sinkSetCache[i]!=NULL) free(sinkSetCache[i]);
    }

    all_sink = NULL;
    all_source = NULL;
    nodeMap = NULL;
    sourceAdjacentArcIndices = NULL;
    sinkAdjacentArcIndices = NULL;
    pdifferenceCache = NULL;

    for(int i=0; i<2; i++)
    {
        arcListCache[i] = NULL;
        nodeListCache[i] = NULL;
        sourceSetCache[i] = NULL;
        sinkSetCache[i] = NULL;
    }
}


static void freeMemoryComplete(void)
/*************************************************************************
freeMemoryComplete
*************************************************************************/
{
	/* destroy breakpoints */
	destroyBreakpoint(firstBreakpoint);
	firstBreakpoint = NULL;

	free(nodeListSuper);
	nodeListSuper = NULL;
	free(arcListSuper);
	arcListSuper = NULL;

    freeCaches();
}

static void freeMemorySolve (void)
{
/*************************************************************************
freeMemorySolve
*************************************************************************/
	uint i;

	for (i=0; i<numNodes; ++i)
	{
		freeRoot (&strongRoots[i]);
	}

	free(strongRoots);
	strongRoots = NULL;

	for (i=0; i<numNodes; ++i)
	{
		if (nodesList[i].outOfTree)
		{
			free(nodesList[i].outOfTree);
			nodesList[i].outOfTree = NULL;
		}
	}

	free(labelCount);
	labelCount = NULL;
}

static void processRoot (Node *strongRoot)
{
/*************************************************************************
processRoot
*************************************************************************/
	Node *temp, *strongNode = strongRoot, *weakNode;
	Arc *out;

	strongRoot->nextScan = strongRoot->childList;

	if ((out = findWeakNode (strongRoot, &weakNode)))
	{
		merge (weakNode, strongNode, out);
		pushExcess (strongRoot);
		return;
	}

	checkChildren (strongRoot);

	while (strongNode)
	{
		while (strongNode->nextScan)
		{
			temp = strongNode->nextScan;
			strongNode->nextScan = strongNode->nextScan->next;
			strongNode = temp;
			strongNode->nextScan = strongNode->childList;

			if ((out = findWeakNode (strongNode, &weakNode)))
			{
				merge (weakNode, strongNode, out);
				pushExcess (strongRoot);
				return;
			}

			checkChildren (strongNode);
		}

		if ((strongNode = strongNode->parent))
		{
			checkChildren (strongNode);
		}
	}

	addToStrongBucket (strongRoot, &strongRoots[strongRoot->label]);
	++ highestStrongLabel;
}


static __inline void quickSort (Arc **arr, const uint first, const uint last)
{
/*************************************************************************
quickSort
*************************************************************************/
	int i=0, j, L=first, R=last, beg[MAX_LEVELS], end[MAX_LEVELS], temp=0;
	Arc *swap;

	if ((R-L) <= 5)
	{// Bubble sort if 5 elements or less
		for (i=R; (i>L); --i)
		{
			swap = NULL;
			for (j=L; j<i; ++j)
			{
				if (isExcess(arr[j]->flow - arr[j+1]->flow) < 0)//(arr[j]->flow < arr[j+1]->flow)
				{
					swap = arr[j];
					arr[j] = arr[j+1];
					arr[j+1] = swap;
				}
			}

			if (!swap)
			{
				return;
			}
		}

		return;
	}

	beg[0]=first;
	end[0]=last+1;

	while (i>=0)
	{
		L=beg[i];
		R=end[i]-1;

		if (L<R)
		{
			swap=arr[L];
			while (L<R)
			{
				while ((isExcess(arr[R]->flow - swap->flow)>=0) && (L<R)) //((arr[R]->flow >= swap->flow) && (L<R))
					R--;

				if (L<R)
				{
					arr[L]=arr[R];
					L++;
				}

				while ((isExcess(arr[L]->flow - swap->flow) <= 0) && (L<R)) //((arr[L]->flow <= swap->flow) && (L<R))
					L++;

				if (L<R)
				{
					arr[R]=arr[L];
					R--;
				}
			}

			arr[L] = swap;

			beg[i+1] = L+1;
			end[i+1] = end[i];
			end[i++] = L;

			if ((end[i]-beg[i]) > (end[i-1]-beg[i-1]))
			{
				temp=beg[i];
				beg[i]=beg[i-1];
				beg[i-1]=temp;
				temp=end[i];
				end[i]=end[i-1];
				end[i-1]=temp;
			}
		}
		else
		{
			i--;
		}
	}
}

static __inline void sort (Node * current)
{
/*************************************************************************
sort
*************************************************************************/
	if (current->numOutOfTree > 1)
	{
		quickSort (current->outOfTree, 0, (current->numOutOfTree-1));
	}
}

static __inline void minisort (Node *current)
{
/*************************************************************************
minisort
*************************************************************************/
	Arc *temp = current->outOfTree[current->nextArc];
	uint i, size = current->numOutOfTree;/*, tempflow = temp->flow;*/
	double tempflow = temp->flow;

	for(i=current->nextArc+1; ((i<size) && (isExcess(tempflow - current->outOfTree[i]->flow) < 0)); ++i)
	{
		current->outOfTree[i-1] = current->outOfTree[i];
	}
	current->outOfTree[i-1] = temp;
}


// static /*ullint*/double checkOptimality (const uint gap)
// {
// /*************************************************************************
// checkOptimality
// *************************************************************************/
// 	uint i, check = 1;
// 	/*ullint*/double mincut = 0;
// 	/*llint*/ double *excess = NULL;
//
// 	Arc *tempArc;
//
// 	excess = (/*llint*/double *) malloc (numNodes * sizeof (double/*llint*/));
// 	if (!excess)
// 	{
// 		printf("Out of memory\n");
// 		exit(0);
// 	}
//
// 	// Pushing depicits from all sink adjacent nodes to the sink
// 	for (i=0; i<nodesList[sink].numOutOfTree; ++i)
// 	{
// 		tempArc = nodesList[sink].outOfTree[i];
// 		if (isExcess(tempArc->from->excess) < 0)
// 		{
// 			if (isExcess((tempArc->from->excess + /*(int)*/ tempArc->flow))  < 0)
// 			{
// 				// Excess is high enough to saturate the arc => Flow on residual arc is zeroed
// 				tempArc->from->excess += /*(int)*/ tempArc->flow;
// 				tempArc->flow = 0;
// 			}
// 			else
// 				// Excess is NOT high enough to saturate the arc => Excess is zeroed
// 			{
// 				tempArc->flow = /*(uint)*/ (tempArc->from->excess + /*(int)*/ tempArc->flow);
// 				tempArc->from->excess = 0;
// 			}
// 		}
// 	}
//
// 	for (i=0; i<numNodes; ++i)
// 	{
// 		excess[i] = 0;
// 	}
//
// 	for (i=0; i<numArcs; ++i)
// 	{
// 		if ((arcList[i].from->label >= gap) && (arcList[i].to->label < gap))
// 		{
// 			mincut += arcList[i].capacity;
// 		}
//
// 		if ((isExcess(arcList[i].flow - arcList[i].capacity)>0) || (isExcess(arcList[i].flow) < 0))
// 		{
// 			check = 0;
// 			printf("Warning - Capacity constraint violated on arc (%d, %d). Flow = %lf, capacity = %lf\n",
// 				   arcList[i].from->number,
// 				   arcList[i].to->number,
// 				   arcList[i].flow,
// 				   arcList[i].capacity);
// 		}
// 		excess[arcList[i].from->number - 1] -= arcList[i].flow;
// 		excess[arcList[i].to->number - 1] += arcList[i].flow;
// 	}
//
// 	for (i=0; i<numNodes; i++)
// 	{
// 		if ((i != (source)) && (i != (sink)))
// 		{
// 			if (isExcess(excess[i]))
// 			{
// 				check = 0;
// 				printf ("Warning - Flow balance constraint violated in node %d. Excess = %lf\n",
// 						i+1,
// 						excess[i]);
// 			}
// 		}
// 	}
//
// 	check = 1;
//
// 	if (isExcess(excess[sink] - mincut) != 0)//(excess[sink] != mincut)
// 	{
// 		check = 0;
// 		printf("Warning - Flow is not optimal - max flow does not equal min cut!\n");
// 	}
//
// // 	if (check)
// // 	{
// // 		printf ("Solution checks as optimal. \t Max Flow: \t %lld\n", mincut);
// // 	}
//
// 	free (excess);
// 	excess = NULL;
// 	return mincut;
// }

// static void decompose (Node *excessNode, const uint source, uint *iteration)
// {
// /*************************************************************************
// decompose
// *************************************************************************/
// 	Node *current = excessNode;
// 	Arc *tempArc;
// 	/*uint*/double bottleneck = excessNode->excess;
//
// 	// Find the bottleneck along a path to the source or on a cycle
// 	for ( ;(current->number != source) && (current->visited < (*iteration)) /*&& (current->nextArc < current->numOutOfTree)*/; // Added by Cheng
// 		 current = tempArc->from)
// 	{
// 		current->visited = (*iteration);
// 		tempArc = current->outOfTree[current->nextArc];
//
// 		if (isExcess(tempArc->flow - bottleneck) < 0) //(tempArc->flow < bottleneck)
// 		{
// 			bottleneck = tempArc->flow;
// 		}
// 	}
//
// 	if (current->number == source) // the DFS reached the source
// 	{
// 		excessNode->excess -= bottleneck;
// 		current = excessNode;
//
// 		// Push the excess all the way to the source
// 		while (current->number != source)
// 		{
// 			tempArc = current->outOfTree[current->nextArc]; // Pick arc going out of node to push excess to
// 			tempArc->flow -= bottleneck; // Push back bottleneck excess on this arc
//
// 			if (isFlow(tempArc->flow)) // If there is still flow on this arc do sort on this current node
// 			{
// 				minisort(current);
// 			} else {
// 				++ current->nextArc;
// 			}
//
// 			current = tempArc->from;
// 		}
// 		return;
// 	}
//
// 	++ (*iteration);
//
// 	// This part is temporarily added by Cheng Lu.
// //	if (current->visited == (*iteration)-1)
// 		bottleneck = current->outOfTree[current->nextArc]->flow;
// /*	else{
// 		current = excessNode;
// 		bottleneck = excessNode->excess;
// 	}*/
//
// 	while (current->visited < (*iteration)) //&& (current->nextArc < current->numOutOfTree)) // Added by Cheng
// 	{
// 		current->visited = (*iteration);
// 		tempArc = current->outOfTree[current->nextArc];
//
// 		if (isExcess(tempArc->flow - bottleneck) < 0)//(tempArc->flow < bottleneck)
// 		{
// 			bottleneck = tempArc->flow;
// 		}
// 		current = tempArc->from;
// 	}
//
// 	++ (*iteration);
//
// 	// This part is temporarily added by Cheng Lu.
// /*	if (current->visited == (*iteration)-1)
// 		;
// 	else{
// 		current = excessNode;
// 	}*/
//
// 	while (current->visited < (*iteration)) //&& (current->nextArc < current->numOutOfTree)) // Added by Cheng
// 	{
// 		current->visited = (*iteration);
//
// 		tempArc = current->outOfTree[current->nextArc];
// 		tempArc->flow -= bottleneck;
//
// 		if (isFlow(tempArc->flow))
// 		{
// 			minisort(current);
// 			current = tempArc->from;
// 		} else {
// 			++ current->nextArc;
// 			current = tempArc->from;
// 		}
// 	}
// }

// static void recoverFlow (const uint gap)
// {
// /*************************************************************************
// recoverFlow
// *************************************************************************/
// 	uint iteration = 1;
// 	uint i, j;
// 	Arc *tempArc;
// 	Node *tempNode;
//
// 	Node **nodePtrArray;
//
// 	if ((nodePtrArray = (Node **) malloc (numNodes * sizeof (Node *))) == NULL)
// 	{
// 		printf("Out of memory\n");
// 		exit(0);
// 	}
//
// 	for (i=0; i < numNodes ; i++)
// 		nodePtrArray[i] = &nodesList[i];
//
// 	// Adding arcs FROM the Source to Source adjacent nodes.
// 	for (i=0; i<nodesList[source].numOutOfTree; ++i)
// 	{
// 		tempArc = nodesList[source].outOfTree[i];
// 		addOutOfTreeNode (tempArc->to, tempArc);
// 	}
//
// 	// Zeroing excess on source and sink nodes
// 	nodesList[source].excess = 0;
// 	nodesList[sink].excess = 0;
//
// 	for (i=0; i<numNodes; ++i)
// 	{
// 		tempNode = &nodesList[i];
//
// 		if ((i == (source)) || (i == (sink)))
// 		{
// 			continue;
// 		}
//
// 		if (tempNode->label >= gap) //tempNode is in SINK set
// 		{
// 			tempNode->nextArc = 0;
// 			if ((tempNode->parent) && (isFlow(tempNode->arcToParent->flow)))
// 			{
// 				addOutOfTreeNode (tempNode->arcToParent->to, tempNode->arcToParent);
// 			}
//
// 			for (j=0; j<tempNode->numOutOfTree; ++j)
// 			{ // go over all sink-set-node's arcs and look for arc with NO flow
// 				if (!isFlow(tempNode->outOfTree[j]->flow))
// 				{	// Remove arc with no flow
// 					-- tempNode->numOutOfTree;
// 					tempNode->outOfTree[j] = tempNode->outOfTree[tempNode->numOutOfTree];
// 					-- j;
// 				}
// 			}
//
// 			sort(tempNode);
// 		}
// 	}
//
// 	for (i=lowestPositiveExcessNode ; i < numNodes ; ++i)
// 	{
// 		tempNode = nodePtrArray[i];
// 		while (isExcess(tempNode->excess) > 0)
// 		{
// 			++ iteration;
// 			decompose(tempNode, source, &iteration);
// 		}
// 	}
//
// 	free(nodePtrArray);
// 	nodePtrArray = NULL;
// }

static void readGraphSuper(double * arcMatrix)
/*************************************************************************
readData
*************************************************************************/
{
	if ((nodeListSuper = (Node *)malloc(numNodesSuper * sizeof(Node))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}

	if ((nodeBreakpoints = (double *)malloc(numNodesSuper * sizeof(double))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}


	if ((arcListSuper = (Arc *)malloc(numArcsSuper * sizeof(Arc))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}

	/* Initialization */
	for (int i = 0; i < numNodesSuper; ++i)
	{
		initializeNode(&nodeListSuper[i], i);
		nodeListSuper[i].originalIndex = i;
        nodeBreakpoints[i] = LAMBDA_HIGH;
	}

    nodeBreakpoints[source] = LAMBDA_LOW;

	for (int i = 0; i < numArcsSuper; ++i)
	{
		initializeArc(&arcListSuper[i]);
	}

	if (LAMBDA_LOW == LAMBDA_HIGH)
	{
		useParametricCut = 0;
	}

	for (int i=0; i < numArcsSuper; ++i)
	{
		int from = (int) arcMatrix[i * 4 + 0];
		int to = (int) arcMatrix[i * 4 + 1];
		double constantCapacity = arcMatrix[ i * 4 + 2 ];
		double multiplierCapacity = arcMatrix[ i * 4 + 3 ];

		arcListSuper[i].constant = constantCapacity;
		arcListSuper[i].multiplier = multiplierCapacity;
		arcListSuper[i].from = &nodeListSuper[from];
		arcListSuper[i].to = &nodeListSuper[to];


        double valueLow = constantCapacity + LAMBDA_LOW*multiplierCapacity;
        double valueHigh = constantCapacity + LAMBDA_HIGH*multiplierCapacity;
        if(valueLow < 0.0 || valueHigh < 0.0){
            roundNegativeCapacity = 1;
        }

		++nodeListSuper[from].numAdjacent;
		++nodeListSuper[to].numAdjacent;

        if ( multiplierCapacity != 0.0 )
        {
            ++ numParametricArcs;
        }

	}


}

static int cmp_arc_breakpoint( const void *l, const void *r)
{
    const Arc *la = *(const Arc**) l;
    const Arc *ra = *(const Arc**) r;

    double a_bp = -(la->constant) / la->multiplier;
    double b_bp = -(ra->constant) / ra->multiplier;


    /*
    printf("c comparing %d->%d , %d->%d, %lf / %lf = %lf and %lf / %lf = %lf\n",
            la->from->number, la->to->number,
            ra->from->number, ra->to->number,
            -(la->constant), la->multiplier, a_bp,
            -(ra->constant), ra->multiplier, b_bp);
    if ( a_bp < b_bp ) return -1;
    if ( a_bp < b_bp ) return 1;
    return 0;
*/
    return (a_bp > b_bp) - (a_bp < b_bp);
}


// Sort parametric arcs list according to their breakpoint ( change in piecewise linear
// function)
static void sortParamArcList()
{
    qsort( paramArcList, numParametricArcs, sizeof(Arc*), cmp_arc_breakpoint);

}

static void debugPrintParamArcList()
{
    printf("c Printing parametric Arcs\n");
    int i;
    for ( i=0; i < numParametricArcs; ++i)
    {
        printf("c %d -b/a = %lf/%lf = %lf\n",
                i,
                -paramArcList[i]->constant,
                paramArcList[i]->multiplier,
                -paramArcList[i]->constant/ paramArcList[i]->multiplier);
    }

}

static void allocateParamArcList()
{
    int i = 0;
    int idx = 0;
	if ((paramArcList= (Arc **)malloc(numParametricArcs * sizeof(Arc*))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}


    for ( i = 0; i < numArcsSuper; ++i)
    {
        uint from = arcListSuper[i].from->number;
        uint to = arcListSuper[i].to->number;

        if ( (from == sourceSuper || to == sinkSuper)
                && arcListSuper[i].multiplier != 0.0)
        {
#ifdef VERBOSE
            printf("c add arc %d->%d with a = %lf b=%lf to parametric arcs\n",
                    from,
                    to,
                    arcListSuper[i].multiplier,
                    arcListSuper[i].constant);
#endif
            paramArcList[idx++] = &(arcListSuper[i]);
        }
    }

}


static void pseudoflowPhase1 (void)
{
/*************************************************************************
pseudoflowPhase1
*************************************************************************/
	Node *strongRoot;
	while ((strongRoot = getHighestStrongRoot ()))
	{
		processRoot (strongRoot);
	}
}

static void prepareOutput (int * numBreakpoints, double ** cuts, double ** breakpoints, int stats[5] )
{
/*************************************************************************
printOutput
*************************************************************************/
	Breakpoint *currentBreakpoint;
	int i;
	int j;

	stats[0] = numArcScans;
	stats[1] = numMergers;
	stats[2] = numPushes;
	stats[3] = numRelabels;
	stats[4] = numGaps;

	/* count num breakpoints */
	*numBreakpoints = 0;
	currentBreakpoint = firstBreakpoint;
	while (currentBreakpoint != NULL)
	{
		++*numBreakpoints;
		currentBreakpoint = currentBreakpoint->next;
	}

	double* breakpointsPointer;
	/* print lambda values */
	if ((breakpointsPointer = (double *)malloc(*numBreakpoints * sizeof(double))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}

	currentBreakpoint = firstBreakpoint;
	for (i = 0; i < *numBreakpoints; i++)
	{
		breakpointsPointer[i] = (double) currentBreakpoint->lambdaValue;
		currentBreakpoint = currentBreakpoint->next;
	}

	*breakpoints = breakpointsPointer;

	/* print values nodes*/
	double* cutsPointer;
	if ((cutsPointer = (double *)malloc( numNodesSuper * sizeof(double))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}

    for (j = 0; j < numNodesSuper; j++)
    {
        cutsPointer[j] = nodeBreakpoints[j];
    }

	*cuts = cutsPointer;
}



static void copyArcNew(CutProblem *problem, int *nodeMap, Arc *old, Arc *new,
	double lambda)
/*************************************************************************
copyArcNew - copy basic info arc and point to new nodes
*************************************************************************/
{
	uint newIndexFrom;
	uint newIndexTo;

	initializeArc(new);
	new->capacity = old->multiplier * lambda + old->constant;

    if (new->capacity < 0)
    {
        if (roundNegativeCapacity)
        {
            new->capacity = 0;
        }
        else
        {
            printf("Negative capacity for lambda equal to %lf. Set roundNegativeCapacity to 1 if the value should be rounded to 0.\n",problem->lambdaValue);
            exit(0);
        }
    }


	/* set start and end node */
	newIndexFrom = nodeMap[old->from->number];
	newIndexTo = nodeMap[old->to->number];
	new->from = &problem->nodeList[newIndexFrom];
	new->to = &problem->nodeList[newIndexTo];

	/* update degree nodes*/
	++ new->from->numAdjacent;
	++ new->to->numAdjacent;
}

static void copyArcAdd(Arc *old, Arc *new, double lambda)
/*************************************************************************
copyArcAdd - update arc by adding another
*************************************************************************/

{
    double additional_capacity = old->multiplier * lambda + old->constant;

    if (additional_capacity < 0)
    {
        if (roundNegativeCapacity)
        {
            additional_capacity = 0;
        }
        else
        {
            printf("Negative capacity for lambda equal to %lf. Set roundNegativeCapacity to 1 if the value should be rounded to 0.\n", lambda);
            exit(0);
        }
    }
    new->capacity +=  additional_capacity;
}

static void destroyProblem(CutProblem *problem, int destroySourceSetIndicator)
/*************************************************************************
destroyProblem - Destruct function for CutProblem struct
*************************************************************************/
{
	problem->sourceSet = NULL;
	problem->sinkSet = NULL;
	problem->nodeList = NULL;
	problem->arcList = NULL;
    if (destroySourceSetIndicator)
    {
        free(problem->optimalSourceSetIndicator);
        problem->optimalSourceSetIndicator = NULL;
    }
}

static void initializeContractedProblem(CutProblem *problem, Node *nodeListProblem,
	uint numNodesProblem, Arc *arcListProblem, uint numArcsProblem,
	const double lambdaValue, char *solutionLow, char *solutionHigh,
    int cacheId)
/*************************************************************************
initializeContractedProblem - Setup problems for parametric cut
*************************************************************************/
{
	uint i, newIndexTo, newIndexFrom;
	uint currentNode = 2;
    uint currentSourceSet = 0;
    uint currentSinkSet = 0;
    uint currentArc = 0;

	/* set cut parameters */
	problem->cutValue = 0;

	/* set solved indicator */
	problem->solved = 0;

	/* initialize optimal cut */
	problem->optimalSourceSetIndicator = NULL;

	/* initialize new lambda value */
	problem->lambdaValue = lambdaValue;
	/* set size of node sets */
	problem->numSourceSet = 0;
	problem->numSinkSet = 0;
	problem->numNodesInList = 2;

	/* allocateSpace for nodeMap */
	if (nodeMap==NULL && (nodeMap = (int *)malloc(numNodesProblem* sizeof(int))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}

    for (i = 0; i < numNodesProblem; i++)
	{
        if (i == sourceSuper || solutionLow[i] == 1)
        {   // Source set nodes
			nodeMap[i] = 0;
            problem->numSourceSet++;
        }
        else if (i == sinkSuper || solutionHigh[i] == 0)
        {
            // sink set nodes
			nodeMap[i] = 1;
            problem->numSinkSet++;
        }
        else
        {
            nodeMap[i] = currentNode;
			++currentNode;
            problem->numNodesInList++;
        }
	}

	/* allocate space for the node sets*/
    int initStructures = nodeListCache[cacheId]==NULL;
	if (nodeListCache[cacheId]==NULL && (nodeListCache[cacheId] = (Node *)malloc(numNodesProblem* sizeof(Node))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}
    problem->nodeList = nodeListCache[cacheId];
	if (sourceSetCache[cacheId]==NULL && (sourceSetCache[cacheId] = (Node *)malloc(numNodesProblem* sizeof(Node))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}
    problem->sourceSet = sourceSetCache[cacheId];
	if (sinkSetCache[cacheId]==NULL && (sinkSetCache[cacheId] = (Node *)malloc(numNodesProblem * sizeof(Node))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}
    problem->sinkSet = sinkSetCache[cacheId];

    // initialize nodes
    for (i = 0; i < problem->numNodesInList; i++)
    {
        initializeNode(&problem->nodeList[i], i);
    }
    for (i = 0; i < problem->numSourceSet; i++)
    {
        initializeNode(&problem->sourceSet[i], i);
    }
    for (i = 0; i < problem->numSinkSet; i++)
    {
        initializeNode(&problem->sinkSet[i], i);
    }
     /* source is always first node */
    problem->nodeList[0].originalIndex = -1; /* indicate artificial source node */
    problem->nodeList[1].originalIndex = -2; /* indicate artificial sink node */

	/* create new node sets*/
	for (i = 0; i < numNodesProblem; i++)
	{
        if(nodeMap[i] > 1){
			problem->nodeList[nodeMap[i]].originalIndex = nodeListProblem[i].originalIndex;
        }else if (nodeMap[i] == 0)
		{
			problem->sourceSet[currentSourceSet].originalIndex = nodeListProblem[i].originalIndex;
			currentSourceSet++;
		}
		else
		{
			problem->sinkSet[currentSinkSet].originalIndex = nodeListProblem[i].originalIndex;
			currentSinkSet++;
		}
	}


    /* allocate space for source and sink arc indices */
	if (sourceAdjacentArcIndices==NULL && (sourceAdjacentArcIndices = (int *)malloc( numNodesProblem*  sizeof(int))) == NULL)

	{
		printf("Out of memory\n");
		exit(0);
	}
	if (sinkAdjacentArcIndices==NULL && (sinkAdjacentArcIndices = (int *)malloc(numNodesProblem  *sizeof(int))) == NULL )
	{
		printf("Out of memory\n");
		exit(0);
	}

	/* initialize indices */
	for (i = 0; i < problem->numNodesInList; i++)
	{
		sourceAdjacentArcIndices[i] = -1;
		sinkAdjacentArcIndices[i] = -1;
	}


	/* allocate space for arcs */
    // Note :  allocating for worst case!
	if (arcListCache[cacheId]==NULL &&(arcListCache[cacheId] = (Arc *)malloc(numArcsProblem * sizeof(Arc))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}
    problem->arcList = arcListCache[cacheId];

    // Allows avoiding mapping two times
    // (for better cache usage)
	/* determine new number of arcs */
	for (i = 0; i < numArcsProblem; i++)
	{
		newIndexFrom = nodeMap[arcListProblem[i].from->number];
		newIndexTo = nodeMap[arcListProblem[i].to->number];

		if (newIndexFrom == newIndexTo || newIndexTo==0 || newIndexFrom==1 || (newIndexFrom == 0 && newIndexTo == 1))
		{
		}
		else if (newIndexFrom == 0)
		{
			if (sourceAdjacentArcIndices[newIndexTo] == -1)
			{
				sourceAdjacentArcIndices[newIndexTo] = currentArc;
				copyArcNew(problem, nodeMap, &arcListProblem[i], &problem->arcList[currentArc], lambdaValue);
				++currentArc;
			}
			else
			{
				copyArcAdd(&arcListProblem[i], &problem->arcList[sourceAdjacentArcIndices[newIndexTo]], lambdaValue);
			}
		}
		else if (newIndexTo == 1)
		{
			if (sinkAdjacentArcIndices[newIndexFrom] == -1)
			{
				sinkAdjacentArcIndices[newIndexFrom] = currentArc;
				copyArcNew(problem, nodeMap, &arcListProblem[i], &problem->arcList[currentArc], lambdaValue);
				++currentArc;
			}
			else
			{
				copyArcAdd(&arcListProblem[i], &problem->arcList[sinkAdjacentArcIndices[newIndexFrom]], lambdaValue);
			}
		}
		else
		{
			copyArcNew(problem, nodeMap, &arcListProblem[i], &problem->arcList[currentArc], lambdaValue);
			++currentArc;
		}
	}

	/* set number of arcs */
	problem->numArcs = currentArc;



}

static void initializeParametricCut(CutProblem *lowProblem, CutProblem *highProblem)
/*************************************************************************
initializeParametricCut - Set up data structures for parametric cut
*************************************************************************/
{
	// disable contraction by passing dummy low/high problem solutions.
    if (all_sink==NULL && (all_sink = (char *)malloc(numNodesSuper *  sizeof(char))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}
	if (all_source==NULL && (all_source = (char *)malloc(numNodesSuper  *sizeof(char))) == NULL )
	{
		printf("Out of memory\n");
		exit(0);
	}
  for (uint i = 0; i < numNodesSuper; i++)
  {
      all_sink[i] = 0;
      all_source[i] = 1;
  }

  /* initialize problem for LAMBDA_LOW */
  initializeContractedProblem(lowProblem, nodeListSuper, numNodesSuper,
		arcListSuper, numArcsSuper,LAMBDA_LOW, all_sink, all_source, 0);

	if (useParametricCut == 1)
	{
		/* initialize problem for LAMBDA_HIGH */
		initializeContractedProblem(highProblem, nodeListSuper, numNodesSuper,
			arcListSuper, numArcsSuper,LAMBDA_HIGH, all_sink, all_source, 1);
	}

    free(all_sink);
    all_sink = NULL;

    free(all_source);
    all_source = NULL;
}

static void addBreakpoint(double lambdaValue, char *sourceSetIndicator)
/*************************************************************************
addBreakpoint - Adds a breakpoint to the linkedlist
*************************************************************************/
{
	Breakpoint *newBreakpoint;
	uint i;


	/* allocate memory for breakpoint*/
	if ((newBreakpoint= (Breakpoint*)malloc(sizeof(Breakpoint))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}

	time( &rawtime );
	info = localtime( &rawtime );
	/* assign values */
	newBreakpoint->lambdaValue = lambdaValue;
	newBreakpoint->next = NULL;

	/* assign space for cut */
    /*
	if ((newBreakpoint->sourceSetIndicator = (uint*)malloc(numNodesSuper * sizeof(uint))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}
    */

	/* copy cut */
	for (i = 0; i < numNodesSuper; i++)
	{
		//newBreakpoint->sourceSetIndicator[i] = sourceSetIndicator[i];
        if ( sourceSetIndicator[i] && nodeBreakpoints[i] > lambdaValue )
        {
            nodeBreakpoints[i] = lambdaValue;
        }
	}

	/* add breakpoint to linkedlist */
	if (lastBreakpoint == NULL)
	{
		/* initialize list */
		firstBreakpoint = newBreakpoint;
		lastBreakpoint = newBreakpoint;
	}
	else
	{
		/* add new element to linked list*/
		lastBreakpoint->next = newBreakpoint;
		/* update head */
		lastBreakpoint = newBreakpoint;
	}
}/*addBreakpoint*/

static void createMemoryStructures( )
/*************************************************************************
createMemoryStructures - creates memory structures
*************************************************************************/
{
	uint from;
	uint to;
	uint i;
	double capacity;

	/* create memory structures */
	for (i=0; i<numNodes; ++i)
	{
		createOutOfTree(&nodesList[i]);
	}

	for (i=0; i<numArcs; i++)
	{
		to = arcList[i].to->number;
		from = arcList[i].from->number;
		capacity = arcList[i].capacity;

		if (!((source == to) || (sink == from) || (from == to)))
		{
			if ((source == from) && (to == sink))
			{
				arcList[i].flow = capacity;
			} else if (to == sink) {
				addOutOfTreeNode(&nodesList[to], &arcList[i]);
			} else {
				addOutOfTreeNode(&nodesList[from], &arcList[i]);
			}
		}
	}

	/* allocate memory for root and label count */
	if ((strongRoots = (Root *)malloc(numNodes * sizeof(Root))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}
	if ((labelCount = (uint *)malloc(numNodes * sizeof(uint))) == NULL)
	{
		printf("Could not allocate memory.\n");
		exit(0);
	}

	/* Initialization of root & labelcount */
	for (i = 0; i<numNodes; ++i)
	{
		initializeRoot(&strongRoots[i]);
		labelCount[i] = 0;
	}
}

static void evaluateCut(CutProblem *problem)
/*************************************************************************
evaluateCut - Evaluates optimal cut parameters for a given problem
*************************************************************************/
{
	uint i;
	int originalIndexFrom;
	int originalIndexTo;
	for (i = 0; i < problem->numArcs; ++i)
	{
		originalIndexFrom = problem->arcList[i].from->originalIndex;
		originalIndexTo = problem->arcList[i].to->originalIndex;
		if ((originalIndexFrom == -1 || problem->optimalSourceSetIndicator[originalIndexFrom] == 1) && (originalIndexTo == -2 || problem->optimalSourceSetIndicator[originalIndexTo] == 0 ) )
		{
		  problem->cutValue += problem->arcList[i].capacity;
		}
	}
}

static void solveProblem(CutProblem *problem, uint maximalSourceSet)
/*************************************************************************
solveProblem - solves a single instance of cut problem
*************************************************************************/
{
	uint i;
	char *tempSourceSet;
	uint nodeCount;

	nodesList = problem->nodeList;
	numNodes = problem->numNodesInList;
	numArcs = problem->numArcs;
	problem->cutValue = 0.0;

    // reset some globals
    highestStrongLabel = 1;
	lowestPositiveExcessNode = 0;


	// handle empty problems
	if (numNodes == 2)
	{
		/* assign nodes to source / sink set */
		if ((problem->optimalSourceSetIndicator = (char *)malloc(numNodesSuper * sizeof(char))) == NULL)
		{
			printf("Out of memory\n");
			exit(0);
		}

		for (i = 0; i < problem->numSourceSet; i++)
		{
			problem->optimalSourceSetIndicator[problem->sourceSet[i].originalIndex] = 1;
		}

		for (i = 0; i < problem->numSinkSet; i++)
		{
			problem->optimalSourceSetIndicator[problem->sinkSet[i].originalIndex] = 0;
		}

		/* determine cut value */
		for (i = 0; i < problem->numArcs; i++)
		{
			if (problem->arcList[i].from->originalIndex == -1 && problem->arcList[i].to->originalIndex == -2)
			{
				problem->cutValue += problem->arcList[i].capacity;
			}
		}

        problem->solved =1;
		return;
	}


	if (maximalSourceSet == 1)
	{
		source = 1;
		sink = 0;

		/* allocate space for reversed arcs */
		if ((arcList = (Arc *)malloc(numArcs * sizeof(Arc))) == NULL)
		{
			printf("Out of memory\n");
			exit(0);
		}

		/* copy arcs such that arcs can be reversed */
		for (i = 0; i < numArcs; i++)
		{
			/* initialize new arc*/
			initializeArc(&arcList[i]);

			// reverse direction
			arcList[i].from = problem->arcList[i].to;
			arcList[i].to = problem->arcList[i].from;

			// assign capacity
			arcList[i].capacity = problem->arcList[i].capacity;
		}
	}
	else
	{
		source = 0;
		sink = 1;

		arcList = problem->arcList;
	}

	// solve
	createMemoryStructures();
	simpleInitialization();
	pseudoflowPhase1();

	/* allocate memory for source set (possibly reversed) */
	nodeCount = problem->numNodesInList + problem->numSourceSet + problem->numSinkSet - 2;
	if ((tempSourceSet = (char *)malloc(nodeCount * sizeof(char))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}

	// retrieve optimal sourceSet for nodes in graph
	if (maximalSourceSet == 1) // reverse assignment to source and sink set
	{
		for (i = 2; i<numNodes; ++i) // start from 2 to ignore artificial source and sink
		{
			if (nodesList[i].label >= numNodes)
			{
				tempSourceSet[nodesList[i].originalIndex] = 0;
            }
			else
			{
				tempSourceSet[nodesList[i].originalIndex] = 1;
			}
		}
	}
	else
	{
		for (i = 2; i<numNodes; ++i) // start from 2 to ignore artificial source and sink
		{
			if (nodesList[i].label >= numNodes)
			{
				tempSourceSet[nodesList[i].originalIndex] = 1;
			}
			else
			{
				tempSourceSet[nodesList[i].originalIndex] = 0;
			}
		}
	}

	// process cut for source set nodes
	for (i = 0; i < problem->numSourceSet; i++)
	{
		tempSourceSet[problem->sourceSet[i].originalIndex] = 1;
	}
	// process cut for sink set nodes
	for (i = 0; i < problem->numSinkSet; i++)
	{
		tempSourceSet[problem->sinkSet[i].originalIndex] = 0;
	}

	// assign cut
	problem->optimalSourceSetIndicator = tempSourceSet;
	evaluateCut(problem);

	if (maximalSourceSet == 1)
	{
		// free if new memory has been allocated for arclist. Memory should not be freed if arcList is taken from the problem
		free(arcList);
		arcList = NULL;
	}

    problem->solved =1;

	//printCutProblem(problem);
	// printf("lambda:%.12lf\n" ,problem->lambdaValue);

	freeMemorySolve();
}

static void differenceSourceSets(char **ppdifference,
	char *lowOptimalSourceIndicator, char *highOptimalSourceIndicator)
{
    if (pdifferenceCache==NULL && (pdifferenceCache = (char *)malloc(numNodesSuper * sizeof(char))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}
    *ppdifference = pdifferenceCache;
    char *pdifference = *ppdifference;
    for (int i = 0; i < numNodesSuper; i++)
    {
      pdifference[i] = highOptimalSourceIndicator[i] - lowOptimalSourceIndicator[i];
    }
}
static void computeCutLinearFunction(char *optimalSourceSetIndicator, double *multiplier, double *constant)
{
    int i;
    int from;
    int to;
    double arc_constant = 0.0;
    double arc_multiplier = 0.0;
    double arc_capacity = 0.0;

    *multiplier = 0.0;
    *constant = 0.0;

    for ( i=0; i < numArcsSuper; i++)
    {
        from = arcListSuper[i].from->originalIndex;
        to = arcListSuper[i].to->originalIndex;
        arc_constant = arcListSuper[i].constant;
        arc_multiplier = arcListSuper[i].multiplier;

        if (optimalSourceSetIndicator[from] &&
                ! optimalSourceSetIndicator[to] )
        {
            *constant += arc_constant;
            *multiplier += arc_multiplier;

        }

    }

}

// Computes intersection of two cut functions that are non-piecewise linear
static double computeLinearIntersect(char *lpS, char *hpS) {

    double a_l, b_l, a_h, b_h;

    a_l = 0.0;
    b_l = 0.0;
    a_h = 0.0;
    b_h = 0.0;
    int i;
    int from;
    int to;
    double arc_constant = 0.0;
    double arc_multiplier = 0.0;

    for ( i=0; i < numArcsSuper; i++)
    {
        from = arcListSuper[i].from->originalIndex;
        to = arcListSuper[i].to->originalIndex;
        arc_constant = arcListSuper[i].constant;
        arc_multiplier = arcListSuper[i].multiplier;

        if (lpS[from] && !lpS[to] )
        {
            b_l += arc_constant;
            a_l += arc_multiplier;
        }

        if (hpS[from] && !hpS[to] )
        {
            b_h += arc_constant;
            a_h += arc_multiplier;
        }

    }

    return -(b_h-b_l)/(a_h-a_l);

}

static double internalCutCapacity(char *optimalSourceSetIndicator) {
    int from, to;
    double arc_capacity;
    double capacity = 0;
		// printf("OG STARTING internalCutCapacity\n");
    for (int i=0; i < numArcsSuper; i++)
    {
        from = arcListSuper[i].from->originalIndex;
        to = arcListSuper[i].to->originalIndex;
        arc_capacity = arcListSuper[i].constant;
        if (optimalSourceSetIndicator[from] == 1
					&& optimalSourceSetIndicator[to] == 0
					&& from != sourceSuper && to != sinkSuper)
        {
            capacity += arc_capacity;
        }
    }
    return capacity;
}

static double evaluateArcFunction(Arc *a, double lambda)
{
    return a->constant + a->multiplier * lambda;

}


double evaluateCutFunction(char *S, double lambda) {
    double sum = 0.0;
    double arc_val;
    int from, to;

    for (int i = 0; i < numArcsSuper; ++i) {
        from = arcListSuper[i].from->originalIndex;
        to   = arcListSuper[i].to->originalIndex;

        if (S[from] && !S[to]) {
            // Arc is in the cut
            arc_val = arcListSuper[i].constant + arcListSuper[i].multiplier * lambda;
            if (arc_val > 0.0) {
                sum += arc_val;
            }
        }
    }

    return sum;
}

// Computes the intersection of the two piecewise-linear functions
// of the low and high cut, computed at lambda_0 and lambda_1 respectively

// We traverse the linear sub-intervals, tracking which arcs linear functions
// are activated (with non zero slope at that interval), and compute the

//
// intersection in linear time respect to the number of arcs






static double computePiecewiseIntersect(char *lpS, char *hpS, double lambda_0, double lambda_1) {
#ifdef DEBUG
    printf("c computing piecewise intersection in [%lf, %lf]\n",  lambda_0, lambda_1);
#endif
    double res = 1.0 / 0.0; // set as infinite by default

    double a_l = 0, b_l = 0, a_h = 0, b_h = 0;
    double arc_constant, arc_multiplier, arc_bp;
    int from, to, i;

#ifdef DEBUG
    char *activated = calloc(numParametricArcs, sizeof(char));
    char *initialized = calloc(numParametricArcs, sizeof(char));
#endif

    // Initialization step
    for (i = 0; i < numArcsSuper; ++i) {
        from = arcListSuper[i].from->originalIndex;
        to = arcListSuper[i].to->originalIndex;
        arc_constant = arcListSuper[i].constant + 0.0;
        arc_multiplier = arcListSuper[i].multiplier + 0.0;

        arc_bp = -arc_constant/arc_multiplier +0.0;

        int active_now = ( (arc_multiplier > 0.0 && arc_bp <= lambda_0)
                || (arc_multiplier < 0.0 && arc_bp > lambda_0)
                || arc_multiplier == 0.0 );

#ifdef DEBUG
        int j = arcListSuper
        for (int j = 0; j < numParametricArcs; ++j) {
            if (arcListSuper + i == paramArcList[j]) {
                initialized[j] = active_now ? 1 : 0;
                activated[j] = active_now ? 1 : 0;
            }
        }
#endif

        if (!active_now) continue;

        if (lpS[from] && !lpS[to]) {
            a_l += arc_multiplier;
            b_l += arc_constant;
        }
        if (hpS[from] && !hpS[to]) {
            a_h += arc_multiplier;
            b_h += arc_constant;
        }
    }

#ifdef VERBOSE
    printf("c Initial coefficients: a_l=%.4lf, b_l=%.4lf | a_h=%.4lf, b_h=%.4lf\n", a_l, b_l, a_h, b_h);
    printf("c Source set lpS: ");
    for (i = 0; i < numNodesSuper; ++i) printf("%d", lpS[i]);
    printf("\n");
    printf("c Source set hpS: ");
    for (i = 0; i < numNodesSuper; ++i) printf("%d", hpS[i]);
    printf("\n");
#endif

    for (i = 0; i < numParametricArcs; ++i) {
        arc_constant = paramArcList[i]->constant + 0.0;
        arc_multiplier = paramArcList[i]->multiplier + 0.0;
        arc_bp = -arc_constant / arc_multiplier + 0.0;
        if (arc_bp > lambda_0) break;
    }

    double lambda_prev = lambda_0 + 0.0;
    double lambda_next;

    for (; i <= numParametricArcs;) {
        lambda_next = (i < numParametricArcs)
            ? -paramArcList[i]->constant / paramArcList[i]->multiplier + 0.0
            : lambda_1 + 0.0;

        if (lambda_next > lambda_1) lambda_next = lambda_1 + 0.0;

        double denom = a_l - a_h;
        if (fabs(denom) > 0.0) {
            double lambda_star = (b_h - b_l) / denom + 0.0;
#ifdef VERBOSE
            printf("c check intersection in range [%.6lf, %.6lf] is %.6lf\n", lambda_prev, lambda_next, lambda_star);
            double cut_l = evaluateCutFunction(lpS, lambda_star);
            double cut_h = evaluateCutFunction(hpS, lambda_star);
            printf("c low cut at %.6lf = %.6lf\n", lambda_star, cut_l);
            printf("c high cut at %.6lf = %.6lf\n", lambda_star, cut_h);
#endif
            if (lambda_star >= lambda_prev && lambda_star < lambda_next) {
#ifdef VERBOSE
                printf("c valid intersection at lambda = %.12lf\n", lambda_star);
#endif
#ifdef DEBUG
                free(activated);
                free(initialized);
#endif
                res = lambda_star;
                return res;
            }
        }

        if (lambda_next == lambda_1) break;
        if (i == numParametricArcs) break;

        while (i < numParametricArcs) {
            arc_constant = paramArcList[i]->constant + 0.0;
            arc_multiplier = paramArcList[i]->multiplier + 0.0;
            arc_bp = -arc_constant / arc_multiplier + 0.0;

            if (arc_bp > lambda_next ) break;

            from = paramArcList[i]->from->originalIndex;
            to = paramArcList[i]->to->originalIndex;

            int in_lp_cut = lpS[from] && !lpS[to];
            int in_hp_cut = hpS[from] && !hpS[to];

            int activating = (arc_multiplier > 0);
            double sign = activating ? 1 : -1;

#ifdef DEBUG
            int arc_index = i;
            if (!initialized[arc_index]) {
                printf("DEBUG WARNING: Arc %d (%d->%d) was not initialized but is now being updated at lambda = %.6lf\n",
                       arc_index, from, to, lambda_prev);
            }
            if (activating && activated[arc_index]) {
                printf("DEBUG ERROR: Activating already activated arc %d (%d->%d)\n", arc_index, from, to);
            }
            if (!activating && !activated[arc_index]) {
                printf("DEBUG ERROR: Deactivating already deactivated arc %d (%d->%d)\n", arc_index, from, to);
            }
            activated[arc_index] = activating ? 1 : 0;
#endif

            if (in_lp_cut) {
                a_l += arc_multiplier * sign;
                b_l += arc_constant * sign;
#ifdef VERBOSE
                printf("c %s low arc (%d->%d): a_l=%.2lf, b_l=%.2lf\n", activating ? "Activating" : "Deactivating", from, to, a_l, b_l);
#endif
            }
            if (in_hp_cut) {
                a_h += arc_multiplier * sign;
                b_h += arc_constant * sign;
#ifdef VERBOSE
                printf("c %s high arc (%d->%d): a_h=%.2lf, b_h=%.2lf\n", activating ? "Activating" : "Deactivating", from, to, a_h, b_h);
#endif
            }

            ++i;
        }

        lambda_prev = lambda_next + 0.0;
    }

#ifdef DEBUG
    free(activated);
    free(initialized);
#endif
    return res;
}



static double computeIntersect(char *difference, double K12)
{
    double constant = K12;
    double multiplier = 0;

    for (int i = 0; i < numArcsSuper; i++)
    {

        if (arcListSuper[i].from->originalIndex == sourceSuper
					&& difference[arcListSuper[i].to->originalIndex] == 1)
        {
            constant += arcListSuper[i].constant;
            multiplier += arcListSuper[i].multiplier;
        }
        else if (arcListSuper[i].to->originalIndex == sinkSuper
					&& difference[arcListSuper[i].from->originalIndex] == 1
					&& roundNegativeCapacity == 0)
        {
            constant -= arcListSuper[i].constant;
            multiplier -= arcListSuper[i].multiplier;
        }
    }

    return constant / (- multiplier);
}

static void parametricCut(CutProblem *lowProblem, CutProblem *highProblem)
/*************************************************************************
parametricCut - Recursive function that solves the parametric cut problem
*************************************************************************/
{

    // determine difference between source sets of cut.
    char *pdifference_low_high;
    differenceSourceSets(&pdifference_low_high, lowProblem->optimalSourceSetIndicator,
			highProblem->optimalSourceSetIndicator);
    uint num_nodes_different_low_high = sum_array_char(pdifference_low_high,
			numNodesSuper);

    double Klow, Khigh, K12;
	/* find lambda value for which the optimal cut functions(expressed as a function of lambda)
	for the lower bound and upper bound problem intersect. */
	if (num_nodes_different_low_high > 0)
	{
        // find intersection using method outlined in Hochbaum 2003 on inverse spanning-tree.


        // if we are rounding negative capacities, use piecewise-linear
        // intersection, else use simpler linear intersection
        //double oldLambdaIntersect = computeIntersect(pdifference_low_high, K12);

        double lambdaIntersect = roundNegativeCapacity ?
            computePiecewiseIntersect(lowProblem->optimalSourceSetIndicator,
                highProblem->optimalSourceSetIndicator,
                lowProblem->lambdaValue,
                highProblem->lambdaValue) :
                computeLinearIntersect(lowProblem->optimalSourceSetIndicator,
                        highProblem->optimalSourceSetIndicator);



        //printf("c intersection old at %lf\n", oldLambdaIntersect);

#ifdef VERBOSE
        printf("c intersection at %lf\n", lambdaIntersect);
#endif


        //lambdaIntersect = oldLambdaIntersect;
        //double valLowCut = evaluateCutFunction(lowProblem->optimalSourceSetIndicator, lambdaIntersect);
        //double valHighCut = evaluateCutFunction(highProblem->optimalSourceSetIndicator, lambdaIntersect);

        //printf("c low cut at %lf = %lf\n", lambdaIntersect, valLowCut);
        //printf("c high cut at %lf = %lf\n", lambdaIntersect, valHighCut);
        //lambdaIntersect = myIntersect;
        // find minimal and maximal source set at lambdaIntersect.
        // Add/subtract TOL to prevent numerical issues.
        CutProblem minimalIntersect;
        initializeContractedProblem(&minimalIntersect, nodeListSuper, numNodesSuper,
					arcListSuper, numArcsSuper,math_max(lambdaIntersect - TOL, LAMBDA_LOW),
					lowProblem->optimalSourceSetIndicator, highProblem->optimalSourceSetIndicator,
                    0);

        solveProblem(&minimalIntersect, 0);
        destroyProblem(&minimalIntersect, 0);

		    CutProblem maximalIntersect;
        initializeContractedProblem(&maximalIntersect, nodeListSuper, numNodesSuper, arcListSuper,
					numArcsSuper,math_min(lambdaIntersect + TOL, LAMBDA_HIGH), minimalIntersect.optimalSourceSetIndicator,
					highProblem->optimalSourceSetIndicator, 1 );

        solveProblem(&maximalIntersect, 0);
        destroyProblem(&maximalIntersect, 0);

        // check if lambdaIntersect is a breakpoint by comparing min and max source set.
        char *pdifference_min_max_intersect;

        differenceSourceSets(&pdifference_min_max_intersect, minimalIntersect.optimalSourceSetIndicator,
					maximalIntersect.optimalSourceSetIndicator);
        uint num_nodes_different_min_max = sum_array_char(pdifference_min_max_intersect, numNodesSuper);

        if (num_nodes_different_min_max > 0 )
        {
            // Intersection is a breakpoint
            addBreakpoint(lambdaIntersect, minimalIntersect.optimalSourceSetIndicator);

        }
    		parametricCut(lowProblem, &minimalIntersect);

    		/* recurse for higher subinterval */
    		parametricCut(&maximalIntersect, highProblem);


        /* call destructor function */
        destroyProblem(&minimalIntersect, 1);
        destroyProblem(&maximalIntersect, 1);

	}
}

void reset_globals()
{
	TOL = 1E-7;
	numNodes = 0;
	numArcs = 0;
	numNodesSuper = 0;
	numArcsSuper = 0;
	source = 0;
	sink = 0;
	highestStrongLabel = 1;

	numArcScans = 0;
	numPushes = 0;
	numMergers = 0;
	numRelabels = 0;
	numGaps = 0;

	nodesList = NULL;
	strongRoots = NULL;
	labelCount = NULL;
	arcList = NULL;
	nodeListSuper = NULL;
	arcListSuper = NULL;
	lowestPositiveExcessNode = 0;

	lastBreakpoint = NULL;
	firstBreakpoint = NULL;

	useParametricCut = 1;
	roundNegativeCapacity = 0;

	LAMBDA_LOW = 0;
	LAMBDA_HIGH = 0;

    freeCaches();
}

int cmpArc(const void *a, const void *b){
    uint mgcNum = 10;
    uint a_from = (((Arc *)a)->from->number) >> mgcNum ;
    uint a_to = (((Arc *)a)->to->number) >> mgcNum ;
    uint b_from = (((Arc *)b)->from->number) >> mgcNum ;
    uint b_to = (((Arc *)b)->to->number) >> mgcNum ;
    if(a_from > b_from) return 1;
    if(b_from > a_from) return -1;
    if(a_to > b_to ) return 1;
    if(b_to > a_to) return -1;
    return 0;
}

void hpf_solve(int numNodesIn, int numArcsIn, int sourceIn, int sinkIn, double * arcMatrix,
	double lambdaRange[2], int roundNegativeCapacityIn, int * numBreakpoints, double ** cuts,
	double ** breakpoints, int stats[5], double times[3] )
/*************************************************************************
main - Main function
*************************************************************************/
{
	reset_globals();

	double readStart, readEnd, initStart, initEnd, solveStart, solveEnd;

	numArcScans = 0;
	numMergers = 0;
	numPushes = 0;
	numRelabels = 0;
	numGaps = 0;

	// printf("NumNodes: %d\n", numNodesIn);
	// printf("NumArcs: %d\n", numArcsIn);
	// printf("Lambda Range: [%lf, %lf]\n", lambdaRange[0], lambdaRange[1]);
	// printf("Round if negative: %d\n", roundNegativeCapacity);
	// printf("Arc matrix:\n");
	// for (int i = 0; i < numArcsIn; ++i)
	// {
	// 	printf("Row %d: [%.2lf, %.2lf, %.2lf, %.2lf]\n", i, arcMatrix[i * 4 + 0 ], arcMatrix[i * 4 + 1 ], arcMatrix[i * 4 + 2 ], arcMatrix[i * 4 + 3 ]);
	// }

	readStart = clock();
	// readInput
	numNodesSuper = numNodesIn;
	numArcsSuper = numArcsIn;
	sourceSuper = (uint) sourceIn;
	sinkSuper = (uint) sinkIn;
	LAMBDA_LOW = lambdaRange[0];
	LAMBDA_HIGH = lambdaRange[1];
	if (LAMBDA_LOW == LAMBDA_HIGH)
		useParametricCut = 0;
	roundNegativeCapacity = 0;
	readGraphSuper( arcMatrix );
	readEnd = clock();

    //printf("c sorting arcs and initializing par cut\n");
	initStart = clock();
    qsort(arcListSuper, numArcsSuper, sizeof(Arc), cmpArc);
    if ( roundNegativeCapacity )
    {
        printf("c will round negative capacities to zero\n");
        allocateParamArcList();
        sortParamArcList();
#ifdef VERBOSE
        debugPrintParamArcList();
#endif
    }
	CutProblem lowProblem;
	CutProblem highProblem;
	initializeParametricCut(&lowProblem,&highProblem);
	initEnd = clock();

	solveStart = clock();
	if (useParametricCut == 1)
	{
        // solve lower bound problem
        solveProblem(&lowProblem, 0);
        destroyProblem(&lowProblem, 0);

        // solve upper bound problem
        solveProblem(&highProblem, 0);
        destroyProblem(&highProblem, 0); // Change on 5/2/2024

        // find breakpoints + recurse
		    parametricCut(&lowProblem, &highProblem);

        // add upper bound as final breakpoint for last interval.
        addBreakpoint(highProblem.lambdaValue, highProblem.optimalSourceSetIndicator);

		/* deallocate memory */
		destroyProblem(&lowProblem, 1);
		destroyProblem(&highProblem, 1);
	}
	else
	{
		printf("Solving problem with lambda value %.15lf\n", lowProblem.lambdaValue);
		solveProblem(&lowProblem,0);
		/* add solution as breakpoint */
		addBreakpoint(lowProblem.lambdaValue, lowProblem.optimalSourceSetIndicator);
		/* deallocate memory */
		destroyProblem(&lowProblem, 1);
	}
	solveEnd = clock();

	times[0] = (readEnd - readStart)/CLOCKS_PER_SEC;
	times[1] = (initEnd - initStart)/CLOCKS_PER_SEC;
	times[2] = (solveEnd - solveStart)/CLOCKS_PER_SEC;

	/* RECOVER FLOW NEEDS TO BE ADAPTED TO DEAL WITH PARAMETRIC ALGORITHM */
	//	recoverFlow( numNodes );
	//	flow = checkOptimality (numNodes);

	prepareOutput(numBreakpoints, cuts, breakpoints, stats);

	printf("Stats: [%d, %d, %d, %d, %d]\n", stats[0],stats[1],stats[2],stats[3],stats[4]);
	printf("times: [%lf, %lf, %lf]\n", times[0],times[1],times[2]);
	printf("Num breakpoints: %d\n", *numBreakpoints);
	// printf("breakpoints:\n");
	// for (int i = 0; i < *numBreakpoints; ++i)
	// {
	// 	printf("Breakpoint: %lf\n", (*breakpoints)[i]);
	// 	for (int j = 0; j < (int)numNodesSuper; j++)
	// 	{
	// 		printf("Cut indicator: %d\n", (*cuts)[i * (int) numNodesSuper + j ]);
	// 	}
	// }


	freeMemoryComplete ();


}
