/*
 * ********NOTE: Instead of B, supply transpose of B.
 *
 *  name    : comp_hgtm.c (MATLAB MEX file)
 *  purpose : computation of hypergraph transition matrix
 *  author  : Christoph Kawan
 *  date    : April 24, 2006
 *
 * edited : Feb, 2020 MST
 * mex -R2018a comp_hgtm_dtControl.c
 * mex -R2018a COPTIMFLAGS="-O3 -Oy- -DNDEBUG" comp_hgtm_dtControl.c
 */
#include "mex.h"
#include "math.h"
#include "stdio.h"
#include "time.h"

#include "matrix.h"
#define DEBUG 0 // 0 or 1

// const unsigned int MAX_B_DEPTH = 16;          /* Maximal accepted value for input argument b_depth */
const int FALSE = 0;
const int TRUE = 1;
typedef struct {
	  mwIndex dim;                            /* dimension */
	  mwSize size;                           /* length of array rows, cols and entries */
	  mwSignedIndex *rows;                          /* row array */
	  mwSignedIndex *cols;                          /* column array */
	  mxDouble *entries;                       /* entries array */
} SparseMatrix;

void printmst(SparseMatrix R){
      mexPrintf("R=\n");
  for(mwSize i=0;i<R.size;i++){
    mexPrintf("    (%d,%d)=%f",R.rows[i],R.cols[i],R.entries[i]);
  }
  mexPrintf("\n");    
}



typedef struct {
	mxDouble *pr;
	mwIndex *ir;
	mwIndex *jc;
	size_t dim;
	mwSize size;
} MATLAB_sparse;
typedef struct _node node;
struct _node {
  mwIndex i;                                     /* number of contained node */
  mwIndex s;                                     /* number of hypernode if > 0 */   
  node *parent;                              /* pointer to parental node */
  node *first_child;                         /* pointer to first child */
  node *brother;                             /* pointer to younger brother */
};
int dtC_nPartition;                        /* number of dtControl partition elements = number of labels */
mxDouble *dtC_partitions;    // dtCon partition labels = A_labels
mxDouble *box_labels;        /*    */
// int b_depth;
SparseMatrix R;                              /* transition matrix of the hypergraph */
MATLAB_sparse B;                             /* transition matrix of the graph */
mwIndex hypernodenum;
mwSize MAX_MATRIXENTRIES;
int max_label;
size_t boxnumber;
// int depth_diff;
node treeroot;                               /* root node of the tree */ 
time_t T1,T2;



void quicksort(mwIndex buf[], int st, int fin)   
/* standard quick sort algorithm */
{
    #if DEBUG
        mexPrintf("quicksort:  Start\n");
    #endif

	int l = st;
	int r = fin;
	mwIndex tmp;
	mwIndex ref = buf[(l+r)/2];
	do {
		while ((buf[l] < ref) && (l < fin)) l++;
		while ((buf[r] > ref) && (r > st)) r--;
		if (l <= r)	{
			tmp = buf[l];
			buf[l] = buf[r];
			buf[r] = tmp;
			l++;
			r--;
		}
	} while (l <= r);
	if (st < r) quicksort(buf,st,r);
	if (l < fin) quicksort(buf,l,fin);
	return;
}
void quicksort2(mwIndex buf[], mxDouble buf2[], int st, int fin) 
/* quicksorts buf and puts the entries of buf2 into the same order */
{
    #if DEBUG
        mexPrintf("quicksort2:  Start\n");
    #endif

	int l = st;
	int r = fin;
    mwIndex tmp;
	mxDouble tmp2;
	mwIndex ref = buf[(l+r)/2];
	do {
		while ((buf[l] < ref) && (l < fin)) l++;
		while ((buf[r] > ref) && (r > st)) r--;
		if (l <= r)	{
			tmp = buf[l];
			buf[l] = buf[r];
			buf[r] = tmp;
                        
            		tmp2 = buf2[l];
            		buf2[l] = buf2[r];
            		buf2[r] = tmp2;
			l++;
			r--;
		}
	} while (l <= r);
	if (st < r) quicksort2(buf,buf2,st,r);
	if (l < fin) quicksort2(buf,buf2,l,fin);
	return;
}
int is_entry(mwIndex r, mwIndex s)
/* returns if the entry (r,s) is contained in R */
{
    #if DEBUG
        mexPrintf("is_entry:  Start\n");
    #endif

   int i;
   if (R.size == 0) return -1;
   for (i=0; i<=R.size-1; i++) {
      if ((R.rows[i] == r) && (R.cols[i] == s)) return i;
   }
   return -1;
}
int is_outgoingarcfrom(mwIndex N)
/* returns if there are outgoing arcs from hypernode N */
{
    #if DEBUG
        mexPrintf("is_outgoingarcfrom: N=%d\n", N);
    #endif

   mwSize i;                                  
   if (R.size == 0) {
       #if DEBUG
            mexPrintf("is_outgoingarcfrom: returns False\n");
       #endif
       return FALSE;
   }
   for (i=0; i<=R.size-1; i++) {
      if (R.rows[i] == N){
          #if DEBUG
                mexPrintf("is_outgoingarcfrom: returns True\n");
          #endif
          return TRUE;
      }
   }
   #if DEBUG
        mexPrintf("is_outgoingarcfrom: returns False\n");
   #endif
   return FALSE;
}
int remove_voids()
/* searches the matrix R for a zero column and deletes it together with the corresponding row */
{
	#if DEBUG
        mexPrintf("remove_voids:  Start\n");
    #endif

    mwSize j;
    int incoming_arcs;
    mwIndex i;
	if (R.dim <= 0) return FALSE;
	for (i=1; i<=R.dim; i++) {
	    incoming_arcs = FALSE;
	    for (j=0; j<=R.size-1; j++) {
			if ((R.cols[j] == i) && (R.entries[j] >= 1)) {
				incoming_arcs = TRUE;
				goto NEXT;
			}		    
		}
NEXT:   if (!incoming_arcs) { 
			for (j=0; j<=R.size-1; j++)	{
				if ((R.rows[j] == i) || (R.cols[j] == i)) { 
					R.rows[j] = -1;
					R.cols[j] = -1;
					R.entries[j] = 0;
				}
				else {
					if (R.rows[j] > i) R.rows[j]--;
					if (R.cols[j] > i) R.cols[j]--;
				}
			}
			R.dim--;
			return TRUE;
		}
	}
	return FALSE;
}
void initnode(node *thisnd, node *_parent, node *_brother, mwIndex _i, mwIndex _s)
/* initialise node */
{
    #if DEBUG
        mexPrintf("initnode:  Start\n");
    #endif

    thisnd->i = _i;
    thisnd->s = _s;
    thisnd->parent = _parent;
    thisnd->brother = _brother;
    thisnd->first_child = NULL;
    return;
}  
void node_add_brother(node *thisnd, int _i, mwIndex _s, node *(*nd))
/* adds a younger brother to node thisnd. *nd is a pointer to the corresponing node */
{
    #if DEBUG
        mexPrintf("node_add_brother:  Start\n");
    #endif

    node *helpnd = thisnd;
    node *pre_hnd = NULL;
    do {
	    pre_hnd = helpnd;
	    helpnd = helpnd->brother;
    } while (helpnd != NULL);
    pre_hnd->brother = (node *) mxMalloc(sizeof(node));
    if (pre_hnd->brother == NULL) {
        mexPrintf("ERROR (node_add_brother): Memory allocation failed.\n");
	return;
    }
    initnode(pre_hnd->brother,thisnd->parent,NULL,_i,_s);
    *nd = pre_hnd->brother;
    return;
}
void node_add_child(node *thisnd, int _i, mwIndex _s, node *(*nd))
/* adds a new child to node thisnd. *nd is a pointer to the corresponding node */
{
    #if DEBUG
        mexPrintf("node_add_child:  Start\n");
    #endif

	if (thisnd->first_child == NULL) {
        thisnd->first_child = (node *) mxMalloc(sizeof(node));
	   if (thisnd->first_child == NULL) {
		   mexPrintf("ERROR (node_add_child): Memory allocation failed.\n");
		   *nd = NULL;
		   return;
	   }
	   initnode(thisnd->first_child,thisnd,NULL,_i,_s);
	   *nd = thisnd->first_child; 	
	   return;	
	}
	else {
	   node_add_brother(thisnd->first_child,_i,_s,nd);
	   return;
	}
}
void node_search_arcless_hypernode(node *thisnd, node *(*nd))
/* searches the tree with root node thisnd for a hypernode with no outgoing arcs */
/* *nd is a pointer to the corresponding node which marks the found hypernode    */
{
    #if DEBUG
        mexPrintf("node_search_arcless_hypernode: START\n");
    #endif
	node *helpnd;
	if (*nd != NULL) return;
	if (((int)thisnd->s) > 0) { 
        #if DEBUG
            mexPrintf("  Searching...thisnd->s=%d\n",thisnd->s);
        #endif
	   if (!is_outgoingarcfrom(thisnd->s)) {
            *nd = thisnd;
            #if DEBUG
                mexPrintf("node_search_arcless_hypernode: returns hypernode=%d\n",thisnd->s);
            #endif
            return;     
	   }
	}
	if (thisnd->first_child == NULL) {
		*nd = NULL;
        #if DEBUG
            mexPrintf("node_search_arcless_hypernode: returns NULL\n");
        #endif
		return;
	}
	helpnd = thisnd->first_child;
	do {
        *nd = NULL;
  	  	node_search_arcless_hypernode(helpnd,nd);
	  	if (*nd != NULL) return;
       	        helpnd = helpnd->brother;
	} while (helpnd != NULL);
  	*nd = NULL;
	return;
}
void node_find_labelled_outgoing_arcs_upwards(node *thisnd, mwIndex *SET, int *setindex, mxDouble label)
/* searches the tree from node thisnd upwards for outgoing arcs with label 'label' and */
/* stores the reached nodes in the array 'SET'                                         */
{
    #if DEBUG
        mexPrintf("node_find_labelled_outgoing_arcs_upwards:  Start\n");
    #endif

   	int k;
    mwIndex j;
	int ADDNODE;
// 	int outgoing_label = ceil((double)thisnd->i / (double)depth_diff);
    if (thisnd->i == -1) {
        #if DEBUG
            mexPrintf("checking label = %f, this node i = %d, thisnd->s=%d\n",label, ((thisnd->i)), thisnd->s);
        #endif
        return;
    }
    mxDouble outgoing_label = box_labels[thisnd->i-1];  // Assuming i>=1
    #if DEBUG
//         int aaa = box_labels[thisnd->i];
        mexPrintf("checking label = %f, this node i = %d, thisnd->s=%d it's label = %f\n",
                label, ((thisnd->i)), thisnd->s, outgoing_label);
	#endif
	if (outgoing_label == label) {
        
        #if DEBUG
//             mexPrintf("Output label = checking label.\n");
//             if (B.jc[thisnd->i-1] == B.jc[thisnd->i]-1) mexPrintf("For B.jc, i and i-1 are equal\n");
            mexPrintf("B.jc[thisnd->i-1] = %d, B.jc[thisnd->i]-1= %d\n",B.jc[thisnd->i-1],B.jc[thisnd->i]-1);
        #endif
        
		for (j=B.jc[thisnd->i-1]; j<=B.jc[thisnd->i]-1; j++)
		{
			ADDNODE = TRUE;
			for (k=0; k<=setindex[0]-1; k++) if (SET[k] == B.ir[j]+1) ADDNODE = FALSE;
			if (ADDNODE) {
				/*SET[setindex[0]++] = B.ir[j] + 1;*/
				SET[setindex[0]] = B.ir[j] + 1;
				setindex[0]++;
			}
            
            #if DEBUG
                mexPrintf("setindex[0]=%d, SET[%d]=%d\n",setindex[0],setindex[0]-1,SET[setindex[0]-1]);
            #endif
		}
	}
 	if (thisnd->parent != NULL)	{
	     node_find_labelled_outgoing_arcs_upwards(thisnd->parent,SET,setindex,label);
	}
	return;
}
int node_is_path(node *thisnd, int pathsize, mwIndex *PATH, node *(*nd))
/* returns if the path 'PATH' with size 'pathsize' is contained in the tree with root thisnd */
{
    #if DEBUG
        mexPrintf("node_is_path:  Start\n");
    #endif

	node *helpnd;
	if (thisnd->first_child == NULL) {
		*nd = NULL;
		return FALSE;
	}
	helpnd = thisnd->first_child;
	do
	{
	    if (helpnd->i == PATH[0]) {
			if (pathsize == 1) { 
				if (((int)helpnd->s) != 0)	{
					*nd = helpnd;
					return TRUE;
				}
				else {
					*nd = NULL;
					return FALSE;
				}
			}
			else {
				if (node_is_path(helpnd,pathsize-1,&PATH[1],nd)) return TRUE;
			}
		}
		helpnd = helpnd->brother;
	} while (helpnd != NULL);
	return FALSE;
}
int node_add_path(node *thisnd, int pathsize, mwIndex *PATH, mwIndex path_num, node *(*nd))
/* adds the path 'PATH' with length 'pathsize' to the tree if not already contained */
/* *nd is a pointer to the end of the path. returns TRUE or FALSE                   */
{
    #if DEBUG
        mexPrintf("node_add_path:  Start\n");
    #endif

    mwIndex _s;
    node *helpnd = NULL;
    if (thisnd->first_child == NULL) {
		if (pathsize == 1) _s = path_num; else _s = 0;
		node_add_child(thisnd,PATH[0],_s,&helpnd);
		if (_s != 0) {
			*nd = helpnd;
			return TRUE;
		}
		else {
			return node_add_path(thisnd->first_child,pathsize-1,&PATH[1],path_num,nd);
		}
	}
	else {
		helpnd = thisnd->first_child;
		do {
			if (helpnd->i == PATH[0]) {
				if (pathsize == 1) {
					if (helpnd->s == path_num) {
						mexPrintf("ERROR (node_add_path): Path to be inserted already exists.\n");
						return FALSE;
					}
					else {
						helpnd->s = path_num;
						*nd = helpnd;
						return TRUE;
					}
				}
				else {
					if (node_add_path(helpnd,pathsize-1, &PATH[1], path_num, nd)) return TRUE;
				}
			}
			helpnd = helpnd->brother;
		} while (helpnd != NULL);
		if (pathsize == 1) _s = path_num; else _s = 0; 
		node_add_child(thisnd,PATH[0],_s,&helpnd);
		if (_s == 0) {
			return node_add_path(helpnd, pathsize-1, &PATH[1], path_num, nd);
		}
		else {
			*nd = helpnd;
			return TRUE;
		}
	}
}
int reallocate_matrix_memory()
/* reallocates memory for transition matrix R */
{
    #if DEBUG
        mexPrintf("reallocate_matrix_memory:  Start\n");
    #endif

	MAX_MATRIXENTRIES += 1024;
	R.rows = (mwSignedIndex *) mxRealloc(R.rows, sizeof(mwSignedIndex)* MAX_MATRIXENTRIES);
	R.cols = (mwSignedIndex *) mxRealloc(R.cols, sizeof(mwSignedIndex)* MAX_MATRIXENTRIES);
	R.entries = (mxDouble *) mxRealloc(R.entries, sizeof(mxDouble)* MAX_MATRIXENTRIES);
	return ((R.rows != NULL) && (R.cols != NULL) && (R.entries != NULL));
}
int tree_add_hypernode(int setsize, mwIndex *SET, node *H)
/* adds a new hypernode to the tree if not already contained */
/* and an arc from (hyper)node H to the new hypernode        */
{
    #if DEBUG
        mexPrintf("tree_add_hypernode:  Start\n");
    #endif
   	int isentry;
	node *nd = NULL;
	node *nd2 = NULL;
   	quicksort(SET,0,setsize-1);
	if (!node_is_path(&treeroot,setsize,SET,&nd2)) {
    	        nd = NULL;
		if (node_add_path(&treeroot, setsize, SET, hypernodenum+1, &nd)) {
			if (nd == NULL) {
				mexPrintf("ERROR (tree_add_hypernode): node_add_path was TRUE but no node was returned.\n");
				return FALSE;
			}
			hypernodenum++;
            #if DEBUG
                mexPrintf("Added hypernode s=%d, {",nd->s);
                for(int ii2=0;ii2<setsize;ii2++)
                    mexPrintf("%d, ",SET[ii2]);
                mexPrintf("}\n");
            #endif
		}
		else {
			mexPrintf("ERROR (tree_add_hypernode): Adding hypernode failed.\n");
			return FALSE;
		}
        R.rows[R.size] = H->s;
        R.cols[R.size] = hypernodenum;
		R.entries[R.size] = 1;
        #if DEBUG
            mexPrintf("R.size=%d, R.rows[R.size]=%d, R.cols[R.size]=%d\n",R.size,R.rows[R.size],R.cols[R.size]);
        #endif
		R.size++;
		if (R.size == MAX_MATRIXENTRIES) {
			if (!reallocate_matrix_memory()) {
				mexPrintf("ERROR (tree_add_hypernode): Reallocation of memory failed.\n");
				return FALSE;
			}
		}
	}
	else { 
       	isentry = is_entry(H->s,nd2->s);
		if (isentry >= 0) { 
			R.entries[isentry]++;
            #if DEBUG
                mexPrintf("Hypernode Nn={");
                for(int ii=0;ii<setsize;ii++)
                    mexPrintf("%d, ",SET[ii]);
                mexPrintf("} already exists; edge from %d to Nn already exists. Increased the number of edges. R.entries[%d] = %f, isentry=%d\n", H->s, isentry, R.entries[isentry],isentry);
            #endif
		}
		else {
   	        R.rows[R.size] = H->s;
			R.cols[R.size] = nd2->s;
			R.entries[R.size] = 1;
            #if DEBUG
                mexPrintf("The Hypernode Nn={");
                for(int ii=0;ii<setsize;ii++)
                    mexPrintf("%d, ",SET[ii]);
                mexPrintf("} exists. ");
                mexPrintf("Added an edge from hypernode %d to %d in the hypergraph.\n    R.size=%d, R.rows[R.size]=%d, R.cols[R.size]=%d\n",H->s,nd2->s,R.size,R.rows[R.size],R.cols[R.size]);
            #endif
			R.size++;
			if (R.size == MAX_MATRIXENTRIES) {
				if (!reallocate_matrix_memory()) {
					mexPrintf("ERROR (hg_tree_add_hypernode): Reallocation of memory failed.\n");
					return FALSE;
				}
			}
		}  
	}
    #if DEBUG
        mexPrintf("tree_add_hypernode: returns TRUE\n");
    #endif
	return TRUE;
}
int tree_build()
/* the hypergraph algorithm */
{
	mwIndex *setH;
	node *helpnd = NULL;
	node *H;
    mwIndex oldnds, oldndi;
	int i,sizeofsetH;
	int HYPERNODES_REMOVED;	
	setH = (mwIndex *) mxMalloc(boxnumber*sizeof(mwIndex));
        if (setH == NULL) {
		mexPrintf("ERROR (tree_build): Memory allocation failed.\n");
		return FALSE;
	}
	R.rows = (mwSignedIndex *) mxMalloc(MAX_MATRIXENTRIES*sizeof(mwSignedIndex));
	R.cols = (mwSignedIndex *) mxMalloc(MAX_MATRIXENTRIES*sizeof(mwSignedIndex));
	R.entries = (mxDouble *) mxMalloc(MAX_MATRIXENTRIES*sizeof(mxDouble));
        if ((R.rows == NULL) || (R.cols == NULL) || (R.entries == NULL)) {
		mexPrintf("ERROR (tree_build): Memory allocation failed.\n");
		return FALSE;
	}
    node_add_child(&treeroot,1,1,&helpnd);
	hypernodenum = 1;
	
    oldndi = -100;
	oldnds = -100;
    
    #if DEBUG
        int tempCount = 0;
    #endif
        
	while (TRUE) {
        #if DEBUG
            tempCount++;
            mexPrintf("Number of times while loop = %d\n",tempCount);
        #endif
            
		H = NULL;
		node_search_arcless_hypernode(&treeroot,&H);
//         #if DEBUG
//             mexPrintf("Search result: hypernode = %d\n",H->s);
//         #endif
		if (H == NULL) goto PRUNING;
		if ((H->i == oldndi) && (H->s == oldnds)) {
			mexPrintf("ERROR (tree_build): Matrix B corrupt.\n");
            
            #if DEBUG
                 mexPrintf("\nhypernodenum = %d, oldndi = %d, oldnds = %d\n", hypernodenum, oldndi, oldnds);
            #endif
            
			return FALSE;
		}
		oldndi = H->i;
		oldnds = H->s;
		for (i=0; i<max_label; i++) {
			sizeofsetH = 0;
            
            #if DEBUG
                mexPrintf("\nChecking for label= %f.\n",dtC_partitions[i]);
            #endif
                
			node_find_labelled_outgoing_arcs_upwards(H,setH,&sizeofsetH,dtC_partitions[i]);
	  	        if (sizeofsetH != 0) {
				if (!tree_add_hypernode(sizeofsetH,setH,H)) 
                    return FALSE;                     
                }	   
		}
	}
PRUNING:
    #if DEBUG
        mexPrintf("Pruning\n");
    #endif
	do {
		R.dim = hypernodenum;
		HYPERNODES_REMOVED = remove_voids();
		hypernodenum = R.dim;
	} while (HYPERNODES_REMOVED);	
    #if DEBUG
        mexPrintf("Freeing memory of setH\n");
    #endif   
	mxFree(setH);
    #if DEBUG
        mexPrintf("Pruning complete\n");
        printmst(R);
    #endif
	return TRUE;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* gateway function between MATLAB and the C computational routine */
{
	size_t N,M;
	mwIndex *colpt, i;
    mxDouble *valpt;
	mwIndex *irs,*jcs;
	mxDouble *sr;
    mxComplexDouble *si;
// 	mwIndex counter;
    mwSize nzsize, nzmax, counter, j;
	int k;
 
	time(&T1);
	if (nrhs != 4) {
		mexErrMsgTxt("Four input arguments expected.");
	}
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}
    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0])*mxGetM(prhs[0]) != 1) {
                mexErrMsgTxt("First input argument must be a scalar.");
	}
    
// 	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1])*mxGetM(prhs[1]) != 1) {
//                 mexErrMsgTxt("Second input argument must be a scalar.");
// 	}
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1]) != 1) {
                mexErrMsgTxt("Second input argument must be a column vector of 'doubles'.");
	}

	
	if (!mxIsSparse(prhs[2])) {
		mexErrMsgTxt("Third input argument must be a sparse matrix.");
	}	
    
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetN(prhs[3]) != 1) {
                mexErrMsgTxt("Fourth input argument must be a column vector of 'doubles'.");
	}
	
	dtC_nPartition = mxGetScalar(prhs[0]);
// 	b_depth = mxGetScalar(prhs[1]);
//         if ((a_depth < 0) || (b_depth < 0) || (a_depth > b_depth)) {
// 		mexErrMsgTxt("0 <= a_depth <= b_depth expected.");
// 	}
// 	if (b_depth > MAX_B_DEPTH) {
//                 mexErrMsgTxt("b_depth too large.");
// 	}
//         max_label = (int)pow(2,a_depth);
    max_label = dtC_nPartition;
	boxnumber = mxGetM(prhs[1]);
//     box_labels = mxGetPr(prhs[1]);
//     dtC_partitions = mxGetPr(prhs[3]);
    dtC_partitions = mxGetDoubles(prhs[3]);
    box_labels = mxGetDoubles(prhs[1]);
    
    #if DEBUG
        mexPrintf("dtC_nPartition = %d, boxnumber= %d, box_labels[0]= %f, box_labels[1]= %f dtC_partitions[0]= %f, dtC_partitions[1]= %f\n",dtC_nPartition,boxnumber,box_labels[0],box_labels[1],dtC_partitions[0],dtC_partitions[1]);
        mexPrintf("sizeof box_labels= %d\n",mxGetM(prhs[1]));
    #endif
        
// 	depth_diff = (int)(boxnumber / max_label);
  
        nzmax = mxGetNzmax(prhs[2]);
    	N = mxGetN(prhs[2]);
    	M = mxGetM(prhs[2]);
// 	B.pr = mxGetPr(prhs[2]);
         B.pr = mxGetDoubles(prhs[2]);
    	B.ir = mxGetIr(prhs[2]);       
    	B.jc = mxGetJc(prhs[2]);
	if (N != M) {
		mexErrMsgTxt("Sparse matrix must be quadratic.");
	}
    	if (N != boxnumber) {
                mexErrMsgTxt("Dimension of sparse matrix supposed to match the number of elements in the column vector.");
	}
	B.dim = N;
	B.size = nzmax;
   	initnode(&treeroot,NULL,NULL,-1,-1);
	hypernodenum = 0;
	R.dim = 0;
	R.size = 0;
	MAX_MATRIXENTRIES = boxnumber;
	if (MAX_MATRIXENTRIES < 1024) MAX_MATRIXENTRIES = 1024;
	
	if (!tree_build()) return;
    colpt = (mwIndex *) mxMalloc(R.dim*sizeof(mwIndex));
	valpt = (mxDouble *) mxMalloc(R.dim*sizeof(mxDouble));
	
	if ((colpt == NULL) || (valpt == NULL)) {
        	mexErrMsgTxt("Allocation of memory failed.");
	}
	plhs[0] = mxCreateSparse(R.dim,R.dim,R.size,mxREAL);
	
	if (plhs[0] == NULL) {
        	mexErrMsgTxt("Allocation of memory failed.");
	}
  
// 	sr  = mxGetPr(plhs[0]);
// 	si  = mxGetPi(plhs[0]);
    sr  = mxGetDoubles(plhs[0]);
	si  = mxGetComplexDoubles(plhs[0]);
	irs = mxGetIr(plhs[0]);
	jcs = mxGetJc(plhs[0]);
	nzsize = R.size;
  
    	if ((sr == NULL) || (irs == NULL) || (jcs == NULL)) {
        	mexErrMsgTxt("Allocation of memory failed.");
	}
	counter = 0;
	for (i=1; i<=R.dim; i++) {
	   k = 0;
	   for (j=0; j<=R.size-1; j++) {
           if (R.cols[j] == i) {
			   colpt[k] = R.rows[j];
			   valpt[k] = R.entries[j];
			   k++;
           }
		}
		if (k == 0) {
			mexErrMsgTxt("Column in sparse matrix could not be found");
		}
		else {     
			quicksort2(colpt,valpt,0,k-1);
 			jcs[i-1] = counter;
			for (j=0; j<=k-1; j++) {
				if (counter >= nzsize) {
					nzsize++;
					mxSetNzmax(plhs[0], nzsize); 
					mxSetDoubles(plhs[0], mxRealloc(sr, nzsize*sizeof(mxDouble)));
					mxSetIr(plhs[0], mxRealloc(irs, nzsize*sizeof(mwIndex)));
					sr  = mxGetDoubles(plhs[0]);
					irs = mxGetIr(plhs[0]);
				}     
				irs[counter] = colpt[j]-1;
				sr[counter] = valpt[j];
				counter++;
			}
		}
	}
	jcs[R.dim] = counter;
//     mexPrintf("\nNow freeing memory of valpt and colpt\n");
	mxFree(valpt);
	mxFree(colpt);
	time(&T2);
	mexPrintf("Total time : %f\n", difftime(T2,T1));
	mexPrintf("Hypernodes : %i\n", hypernodenum);
	return;
}
