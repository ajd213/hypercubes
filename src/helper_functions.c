/* Helper functions for the module hypercubes. E.g. for operating the stack, for checking
user input, and for computing common mathematical functions, like binomial coefficients. */
#define NO_IMPORT_ARRAY
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>


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
        printf("N must be between 1 and 32.\n");
        return false;
    }

    if (p < 0 || p > 1)
    {
        printf("p is a probability, and must satisfy 0 <= p <= 1.\n");
        return false;
    }

    if (NR < 1)
    {
        printf("Please use one or more realisations: NR > 1.\n");
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
            printf("Error with realloc in stack!\n");
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
