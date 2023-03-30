#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h> 
#include "functions.h"


// push site to the stack
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

// pop site from the stack
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

// unsigned long integer exponentiation
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

// binomial Coefficient C(n, k)
// inefficient, but used infrequently
ul binomialCoeff(ul n, ul k)
{
    // Base Cases
    if (k > n)
        return 0;

    if (k == 0 || k == n)
        return 1;
 
    // Pascal's triangle
    return binomialCoeff(n - 1, k - 1) + binomialCoeff(n - 1, k);
}

// reset the visited array to all be false
void reset_visited(bool visited[], ul length)
{
    for (ul i = 0; i < length; i++)
    {
        visited[i] = false;
    }
}

// test the stack by putting all of the sites in it
void test_stack(stack *s)
{
    for (ul i = 0; i < s->length; i++)
    {
        s->sites[i] = i;
    }
    s->top = s->length - 1;
}

// binary search an ordered list of sites
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

// populate the sites of the XXZ graph
void populate_sites_XXZ(ul *sites, ul N, int UP)
{
    ul full_HS = intpower(2, N);
    for (ul i = 0, counter = 0; i < full_HS; i++)
    {
        if (__builtin_popcount(i) == UP)
        {
            sites[counter] = i;
            counter++;
        }
    }
}
