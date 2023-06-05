/* Helper functions for the module hypercubes. E.g. for operating the stack, for checking
user input, and for computing common mathematical functions, like binomial coefficients. */
#define NO_IMPORT_ARRAY
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <limits.h>

extern gsl_rng *RNG;

/*
 * Function:  Hamming_distance
 * --------------------
 * Compute the Hamming distance between two basis sites (integers). This is 
 * the number of bits different between the two sites, i.e., the minimum number
 * of spin-flips required to go from one to the other. 
 * 
 * args: two positive Python integers
 * 
 * returns: a pointer of type PyObject * to the Hamming distance
 */
PyObject *Hamming_distance(PyObject *self, PyObject *args)
{
    PyObject *state1, *state2;

    // parse and check arguments
    if (!PyArg_ParseTuple(args, "OO", &state1, &state2)) goto error;

    ul s1 = pyobject_to_ul(state1);
    ul s2 = pyobject_to_ul(state2);
    if (PyErr_Occurred()) goto error;

    int hd = __builtin_popcount(s1 ^ s2);

    return PyLong_FromLong(hd);

    error:
        PyErr_SetString(PyExc_RuntimeError, "Fatal error occurred in Hamming_distance");
        return NULL;
}

PyObject *RNG_test(PyObject *self, PyObject *args)
{
    int number = gsl_rng_uniform_int(RNG, INT_MAX);
    return PyLong_FromLong(number);

}

/*
 * Function:  CArrayToNumPyArray
 * --------------------
 *  Wrap a pointer to clusters data in a NumPy array.
 * 
 * arr: a pointer to an array of uls
 * length: length of array
 * 
 * returns: a pointer of type PyObject * to the array
 */
PyObject *CArrayToNumPyArray(ul *arr, ul length)
{
    npy_intp dims[] = {length};
    PyObject *numpy_array = PyArray_SimpleNewFromData(1, dims, NPY_ULONGLONG, (void *)arr);

    if (!numpy_array)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error: Unable to create NumPy array in CArrayToNumPyArray.");
        return NULL;
    }

    return numpy_array;
}

/*
 * Function:  pyobject_to_ul
 * --------------------
 * Take any python integer (int, np.int64, etc.) and parse it as a C unsigned long
 * 
 * positive_int: a pointer to a PyObject representing a positive integer
 * 
 * returns: N, an unsigned long
 */
ul pyobject_to_ul(PyObject *positive_int)
{
    // should be followed by if (PyErr_Occurred()) {...} to catch overflow error
    ul N;

    // Convert the input to a Python int
    PyObject *py_int = PyNumber_Long(positive_int);
    if (!py_int)
    {
        PyErr_SetString(PyExc_TypeError, "Invalid input type. N must be a positive integer.");
    }

    // Convert Python int to C unsigned long
    N = PyLong_AsUnsignedLong(py_int);

    Py_DECREF(py_int);

    return N;
}


/*
 * Function:  check_args
 * --------------------
 *  check user arguments N, NR and p.
 *
 *  N: the dimension of the hypercube: 1 <= N <= 32
 *  NR: the Number of Realisations: number of clusters to grow. NR >= 1
 *  p: the percolation concentration. 0 <= p <= 1
 *
 *  returns: true if arguments are OK, else false.
 */
bool check_args(ul N, ul NR, float p)
{
    if (N < 1 || N > 32)
    {
        PyErr_SetString(PyExc_ValueError, "Invalid input arguments! N must be between 1 and 32.");
        return false;
    }

    if (p < 0 || p > 1)
    {
        PyErr_SetString(PyExc_ValueError, "Invalid input arguments! p is a probability, and must satisfy 0 <= p <= 1.");
        return false;
    }

    if (NR < 1)
    {
        PyErr_SetString(PyExc_ValueError, "Invalid input arguments! Please use one or more realisations: NR > 1.");
        return false;
    }
    return true;
}

/*
 * Function:  setup_stack
 * --------------------
 *  allocate memory for the stack and initialise values.
 *
 *  NH: the Hilbert space dimension / the dimension of the graph / the total number of nodes
 *
 *  returns: a pointer (stack *) to a stack, or NULL if malloc fails.
 */
stack *setup_stack(ul NH)
{
    stack *s = malloc(sizeof(stack));
    if (!s)
    {
        return NULL;
    }
    s->top = 0;
    s->length = NH;
    s->NH = NH;
    s->sites = malloc(s->length*sizeof(ul));
    if (!s->sites)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error setting up stack");
        return NULL;
    }

    return s;
}

/*
 * Function:  push
 * --------------------
 * push a site to the stack, and re-size the stack with realloc() if necessary.
 * Stack struct defined in header file.
 *
 *  s: a pointer to a stack, defined in functions.h
 *  site: the graph site to be pushed to the stack
 *
 *  returns: 0 if push successful, 1 otherwise.
 */
int push(stack *s, ul site)
{
    if (s->top == s->length - 1)
    {
        // Overflow! Increase the length of the stack.
        ul *new_stackmem = realloc(s->sites, (s->length + s->NH)*sizeof(ul));
        if (new_stackmem == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError, "Error with realloc in push().");
            return 1;
        }
        else
        {
            s->sites = new_stackmem;
            s->length += s->NH;
        }
    }
    s->sites[s->top] = site;
    s->top++;

    // pushed item
    return 0;
}

/*
 * Function:  pop
 * --------------------
 * pop a site from the stack. Throw an error if the stack is empty.
 *
 *  s: a pointer to a stack, defined in functions.h
 *  error: a pointer to an error flag
 *
 *  returns: the site popped from the stack, if successful. Otherwise 0 (and error raised).
 */
ul pop(stack *s, int *error)
{
    if (s->top == 0)
    {
        // cannot pop if the stack is empty
        *error = -1;
        return 0;
    }

    s->top--;
    return s->sites[s->top];
}

/*
 * Function:  reset_visited
 * --------------------
 * loop through visited[], which tracks whether each site has been visited by the DFS algorithm.
 *
 *  visited: a list of bools equal to NH, the number of nodes in the graph.
 *  length: the length of visited (= NH)
 */
void reset_visited(bool visited[], ul length)
{
    for (ul i = 0; i < length; i++)
    {
        visited[i] = false;
    }
}

/*
 * Function:  index_site
 * --------------------
 * search an ordered list of sites for a particular site using a binary search.
 *
 *  sites: a pointer to an ordered list of sites
 *  site: the site to search for
 *  left: the left index of the list (use 0 when calling)
 *  right: the right index of the list (use len(sites) - 1 when calling)
 *  idx_flag: a pointer to a flag which is set to -1 if the site cannot be found
 *
 *  returns: the index of the site, if found. If not found, return 0 and set idx_flag to -1.
 */
ul index_site(ul *sites, ul site, ul left, ul right, int *idx_flag)
{
    if (right >= left) {
        int mid = left + (right - left) / 2;
 
        if (sites[mid] == site)
            return mid;
 
        if (sites[mid] > site)
            return index_site(sites, site, left, mid - 1, idx_flag);
 
        return index_site(sites, site, mid + 1, right, idx_flag);
    }
    // not found
    *idx_flag = -1;
    return 0;
}

/*
 * Function:  intpower
 * --------------------
 * recursive integer exponentiation for ul (unsigned long) type. Used seldomly.
 *
 *  base: the number to be exponentiated
 *  exponent: the power which the base is raised to
 *
 *  returns: base ** exponent (base to the power of exponent)
 */
ul intpower(ul base, ul exponent)
{
    if (exponent == 1)
    {
        return base;
    }
    else
    {
        return base * intpower(base, exponent - 1);
    }
}

/*
 * Function:  binomialCoeff
 * --------------------
 * Efficiently compute binomial Coefficient C(n, r). Note that 
 * we must evaluate the k! term `backwards' to retain exact divisibiity.
 * O(r) time complexity and O(1) space complexity.
 *
 *  n: integer in C(n, r)
 *  r: integer in C(n, r)
 *
 *  returns: C(n, r) = n!/(r!(n-r)!)
 */
ul binomialCoeff(ul n, ul r)
{
    ul result = 1;
 
    // Shorten the loop: C(n, r) = C(n, n-r)
    if (r > n - r)
        r = n - r;

    for (ul i = 0; i < r; i++) {
        result *= (n - i); // n, (n-1), ..., (n-r+1)
        result /= (i + 1); // 1, 2, ..., r
    }
 
    return result;
}

/*
* Function: fibonacci
* -----------------------
* Efficiently compute the n-th Fibonacci number, indexed from zero (0,1,1,2,3,...)
*
* n : the index of the Fibonacci number to calculate
* 
* returns: the n-th Fibonacci number
*/
ul fibonacci(ul n)
{
    assert(n >= 0);

    if (n == 0)
    {
        return 0;
    }

    ul a = 0;
    ul b = 1;
    while (n-- > 1) {
        int t = a;
        a = b;
        b += t;
    }
    return b;
}


queue *setup_queue(ul length)
{
    queue *q = malloc(sizeof(queue));
    if (!q) return NULL;

    q->writeIdx = 0;
    q->readIdx = 0;
    q->length = length + 1; // one extra value needed
    q->sites = malloc(q->length*sizeof(ul));

    if (!q->sites) 
    {
        PyErr_SetString(PyExc_RuntimeError, "Error setting up queue!");
        return NULL;
    }

    return q;
}

void enqueue(queue *q, ul item, int *err)
{
    if ((q->writeIdx + 1) % q->length == q->readIdx)
    {
        // buffer is full, avoid overflow
        PyErr_SetString(PyExc_RuntimeError, "Error! Queue full.");
        *err = 1;
        return;
    }
    q->sites[q->writeIdx] = item;
    q->writeIdx = (q->writeIdx + 1) % q->length;
}

ul dequeue(queue *q, int *err)
{
    if (q->readIdx == q->writeIdx)
    {
        // empty
        PyErr_SetString(PyExc_RuntimeError, "Error! Queue empty!");
        *err = 2;
        return 1;
    }
    ul value = q->sites[q->readIdx];
    q->readIdx = (q->readIdx + 1) % q->length;
    return value;
}

bool empty(queue *q)
{
    if (q->readIdx == q->writeIdx) return true;
    return false;
}
