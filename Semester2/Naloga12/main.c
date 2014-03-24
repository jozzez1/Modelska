#include <stdio.h>
#include <assert.h>
#include "navier_stokes.h"

int main (void)
{
    unsigned int N = 200;
    assert (!(200 & 1));

    N++;
    return 0;
}
