/******************************************************************************
 * src/SortAlgo.cpp
 *
 * Implementations of many sorting algorithms.
 *
 * Note that these implementations may not be as good/fast as possible. Some
 * are modified so that the visualization is more instructive.
 *
 * Futhermore, some algorithms are annotated using the mark() and watch()
 * functions from SortArray. These functions add colors to the illustratation
 * and thereby makes the algorithm's visualization easier to explain.
 *
 ******************************************************************************
 * The algorithms in this file are copyrighted by the original authors. All
 * code is freely available.
 *
 * The source code added by Timo Bingmann and all modifications are
 * copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "SortAlgo.h"

#include <algorithm>
#include <numeric>
#include <limits>
// #include <iostream>
#include <inttypes.h>
#include <queue>

typedef ArrayItem value_type;

// inversion count limit for iterator instrumented algorithms
const unsigned int inversion_count_instrumented = 512;

const struct AlgoEntry g_algolist[] =
{
    { _("Selection Sort"), &SelectionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Insertion Sort"), &InsertionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Binary Insertion Sort"), &BinaryInsertionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Merge Sort"), &MergeSort, UINT_MAX, 512,
      _("Merge sort which merges two sorted sequences into a shadow array,"
        "and then copies it back to the shown array.") },
    { _("Merge Sort (iterative)"), &MergeSortIterative, UINT_MAX, 512,
      _("Merge sort variant which iteratively merges "
        "subarrays of sizes of powers of two.") },
    { _("Quick Sort (LR ptrs)"), &QuickSortLR, UINT_MAX, UINT_MAX,
      _("Quick sort variant with left and right pointers.") },
    { _("Quick Sort (LL ptrs)"), &QuickSortLL, UINT_MAX, UINT_MAX,
      _("Quick sort variant from 3rd edition of CLRS: two pointers on left.") },
    { _("Quick Sort (ternary, LR ptrs)"), &QuickSortTernaryLR, UINT_MAX, UINT_MAX,
      _("Ternary-split quick sort variant, adapted from multikey quicksort by "
        "Bentley & Sedgewick: partitions \"=<?>=\" using two pairs of pointers "
        "at left and right, then copied to middle.") },
    { _("Quick Sort (ternary, LL ptrs)"), &QuickSortTernaryLL, UINT_MAX, UINT_MAX,
      _("Ternary-split quick sort variant: partitions \"<>?=\" using two "
        "pointers at left and one at right. Afterwards copies the \"=\" to middle.") },
    { _("Quick Sort (dual pivot)"), &QuickSortDualPivot, UINT_MAX, UINT_MAX,
      _("Dual pivot quick sort variant: partitions \"<1<2?>\" using three pointers, "
        "two at left and one at right.") },
    { _("Bubble Sort"), &BubbleSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Cocktail Shaker Sort"), &CocktailShakerSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Gnome Sort"), &GnomeSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Comb Sort"), &CombSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Shell Sort"), &ShellSort, UINT_MAX, 1024,
      _("Uses Robert Sedgewick gaps.") },
    { _("Heap Sort"), &HeapSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Smooth Sort"), &SmoothSort, UINT_MAX, 1024,
      wxEmptyString },
    { _("Odd-Even Sort"), &OddEvenSort, UINT_MAX, 1024,
      wxEmptyString },
    { _("Recursive Batcher's Bitonic Sort"), &BitonicSortRe, UINT_MAX, UINT_MAX, _("Older version.") },
    { _("Iterative Batcher's Bitonic Sort"), &BitonicSortIt, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Recursive Batcher's Odd-Even Merge Sort"), &BatcherSortNetworkRe, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Iterative Batcher's Odd-Even Merge Sort"), &BatcherSortNetworkIt, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Cycle Sort"), &CycleSort, 512, UINT_MAX,
      wxEmptyString },
    { _("Radix Sort (LSD)"), &RadixSortLSD, UINT_MAX, 512,
      _("Least significant digit radix sort, which copies item into a shadow "
        "array during counting.") },
    { _("Radix Sort (MSD)"), &RadixSortMSD, UINT_MAX, UINT_MAX,
      _("Most significant digit radix sort, which permutes items in-place by walking cycles.") },
    { _("std::sort (gcc)"), &StlSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::stable_sort (gcc)"), &StlStableSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("std::sort_heap (gcc)"), &StlHeapSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Tim Sort"), &TimSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Block Merge Sort (WikiSort)"), &WikiSort, UINT_MAX, inversion_count_instrumented,
      _("An O(1) place O(n log n) time stable merge sort.") },
    { _("Bogo Sort"), &BogoSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Bozo Sort"), &BozoSort, 10, UINT_MAX,
      wxEmptyString },
    { _("Stooge Sort"), &StoogeSort, 256, inversion_count_instrumented,
      wxEmptyString },
    { _("Slow Sort"), &SlowSort, 128, inversion_count_instrumented,
      wxEmptyString },
    { wxEmptyString, &PairwiseSortingIt, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Recursive Pairwise Sorting Network"), &PairwiseSortingRe, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Iterative Pairwise Sorting Network"), &PairwiseSortingIt, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Recursive Bose-Nelson Sorting Network"), &BoseNelsonSortingRe, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Iterative Bose-Nelson Sorting Network"), &BoseNelsonSortingIt, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Unbalanced Tree Sort"), &UnbalancedTreeSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Stable Selection Sort"), &StableSelectionSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Stable Quick Sort"), &StableQuickSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Stable Stooge Sort"), &StableStoogeSort, 256, inversion_count_instrumented,
      wxEmptyString },
    { _("Iterative Quick Sort (LL ptrs)"), &IterativeQuickSortLL, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Iterative Quick Sort (LR ptrs)"), &IterativeQuickSortLR, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { _("Iterative Stooge Sort"), &IterativeStoogeSort, 256, inversion_count_instrumented,
      wxEmptyString },
    { _("Gravity Sort"), &GravitySort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Less Bogo Sort"), &LessBogoSort, 20, UINT_MAX,
      wxEmptyString },
    { _("Reverse Insertion Sort"), &ReverseInsertionSort, UINT_MAX, 199,
      wxEmptyString },
    { _("Bingo Sort"), &BingoSort, UINT_MAX, UINT_MAX,
      _("A variant of selection sort for equal items.") },
    { _("Proxmap Sort"), &ProxmapSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Spaghetti Sort"), &SpaghettiSort, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Flipped Min Heap Sort"), &FlippedMinHeapSort, UINT_MAX, UINT_MAX,
      wxEmptyString },
    { wxEmptyString, &PairwiseSortingIt, UINT_MAX, inversion_count_instrumented,
      wxEmptyString },
    { _("Circle Sort"), &CircleSort, UINT_MAX, UINT_MAX,
      _("Custom sort.") },
    { _("Iterative Circle Sort"), &IterativeCircleSort, UINT_MAX, UINT_MAX,
      _("Custom sort.") },
    { _("Balanced Sorting Network"), &BalancedSortingNetwork, UINT_MAX, UINT_MAX,
      _("Custom sort.") },
    { _("Pancake Sort"), &PancakeSort, UINT_MAX, 256,
      _("Custom sort.") },
    { _("Quasi-Pancake Sort"), &QuasiPancakeSort, UINT_MAX, 256,
      _("Custom sort.") },
    { _("Binary Quasi-Pancake Sort"), &BinaryQuasiPancakeSort, UINT_MAX, 256,
      _("Custom sort.") },
    { _("Shove Sort"), &ShoveSort, 256, UINT_MAX,
      _("Custom sort.") },
//  { _("Improved Weave Merge Sort"), &ImprovedWeaveMergeSort, UINT_MAX, inversion_count_instrumented,
//    _("Custom sort.") },
    { _("Iterative Linked Quick Sort (LL ptrs)"), &IterativeLinkedQuickSortLL, UINT_MAX, UINT_MAX,
      _("Custom sort. Uses a queue to keep track of the split ranges.") },
    { _("Iterative Linked Quick Sort (LR ptrs)"), &IterativeLinkedQuickSortLR, UINT_MAX, UINT_MAX,
      _("Custom sort. Uses a queue to keep track of the split ranges.") },
    { _("More Optimized Bubble Sort"), &SmarterBubbleSort, UINT_MAX, UINT_MAX,
      _("Custom sort by Anonymous0726.") },
    { _("Slope Sort"), &SlopeSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Slide Sort"), &SlideSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Grate Sort"), &GrateSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Reverse Grate Sort"), &ReverseGrateSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Cocktail Grate Sort"), &CocktailGrateSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Quad Stooge Sort"), &QuadStoogeSort, 512, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Balance Sort"), &BalanceSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Linked Triangular Heap Sort"), &TriangularHeapSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Checkerboard Heap Sort"), &CheckerboardHeapSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Hyper Stooge Sort"), &HyperStoogeSort, 8, inversion_count_instrumented,
      _("Custom sort by EilrahcF.") },
    { _("Binary Sort"), &BinarySort, 8, inversion_count_instrumented,
      _("Custom sort by EilrahcF.") },
    { _("Bogo Bogo Bogo Sort"), &BogoBogoBogoSort, 10, UINT_MAX,
      _("Custom sort by EilrahcF. This version uses Bogo instead of Bogobogo instead.") },
    { _("Room Sort"), &RoomSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF. Each block is of size sqrt(n).") },
    { _("Wiggle Sort"), &WiggleSort, UINT_MAX, UINT_MAX,
      _("Custom sort by EilrahcF.") },
    { _("Optimized Stable Quick Sort"), &OptiStableQuickSort, UINT_MAX, UINT_MAX,
      _("Custom sort by aphitorite.") },
    { _("Rotate Radix LSD Sort"), &RotateRadixLSD, 512, UINT_MAX, 
      _("Custom sort by aphitorite.") },
    { _("Rotate Radix MSD Sort"), &RotateRadixMSD, 512, UINT_MAX, 
      _("Custom sort by aphitorite.") },
    { _("Awkward Sort"), &AwkwardSort, 256, UINT_MAX, 
      _("Custom sort by aphitorite.") },
    { _("Odd-Even Transposition Merge Sort"), &OddEvenTransMergeSort, UINT_MAX, UINT_MAX, 
      _("Custom sort by aphitorite.") },
    { _("Bubble Scan Quicksort"), &BubbleScanQuicksort, UINT_MAX, UINT_MAX,
      _("Custom sort by aphitorite and thatsOven.") },
    { _("Radix LSD Stooge Sort"), &RadixLSDStoogeSort, 256, inversion_count_instrumented,
      _("Custom sort by thatsOven.") },
    { _("Modulo Sort"), &ModuloSort, 512, UINT_MAX,
      _("Custom sort by McDude_73.") },
    { _("Jump Down Sort"), &JumpDownSort, UINT_MAX, UINT_MAX,
      _("Custom sort by fungamer2.") },
    { _("Cocktail Shell Sort"), &CocktailShellSort, UINT_MAX, 1024,
      _("Custom sort by fungamer2. Each gap decreases by a factor of sqrt(n).") },
    { _("Ternary Slow Sort"), &TernarySlowSort, 128, inversion_count_instrumented,
      _("Custom sort by fungamer2.") },
    { _("Odd Even Bogo Sort"), &OddEvenBogoSort, 256, UINT_MAX,
      _("Custom sort by fungamer2.") },
    { _("Circloid Sort"), &CircloidSort, UINT_MAX, UINT_MAX,
      _("Custom sort by yuji.") },
    { _("Kiyomi Sort"), &KiyomiSort, UINT_MAX, UINT_MAX,
      _("Custom sort by yuji.") },
    { _("Kinnow Sort"), &KinnowSort, UINT_MAX, UINT_MAX,
      _("Custom sort by yuji.") },
    { _("Chinotto Sort"), &ChinottoSort, UINT_MAX, UINT_MAX,
      _("Custom sort by yuji.") },
    { _("Mandarin Orange Sort"), &MandarinOrangeSort, UINT_MAX, UINT_MAX,
      _("Custom sort by yuji.") },
    { _("Y-Slow Sort"), &YSlowSort, 64, inversion_count_instrumented,
      _("Custom sort by yuji.") },
    { _("Float Sort"), &FloatSort, UINT_MAX, UINT_MAX,
      _("Custom sort by Lance.") },
    { _("Odd Even Sort (Base 3)"), &OddEvenBase3, UINT_MAX, UINT_MAX,
      _("Custom sort by Lance.") },
    { _("Snuffle Sort"), &SnuffleSort, 32, inversion_count_instrumented,
      _("Custom sort by _fluffyy.") },
    { _("Quick Sort (LLL pointers)"), &QuickSortLLL, UINT_MAX, UINT_MAX,
      _("Custom sort by u-ndefined.") },
    { _("Radix MSD Stooge Sort"), &RadixMSDStoogeSort, 256, inversion_count_instrumented,
      _("Custom sort by u-ndefined.") }
};

const size_t g_algolist_size = sizeof(g_algolist) / sizeof(g_algolist[0]);

const struct AlgoEntry* g_algolist_end = g_algolist + g_algolist_size;


// ****************************************************************************
// *** Some useful functions for algorithms

void InsertionSortExtra(SortArray& A, ssize_t lo, ssize_t hi)
{
    for (size_t i = lo + 1; i < hi; ++i)
    {
        ArrayItem tmp, key = A[i];
        A.mark(i);

        ssize_t j = i - 1;
        while (j >= lo && (tmp = A[j]) > key)
        {
            A.set(j + 1, tmp);
            j--;
        }
        A.set(j + 1, key);

        A.unmark(i);
    }
}
 
// ****************************************************************************

bool CheckSorted(SortArray& A)
{
    size_t i;
    ArrayItem prev = A[0];
    A.mark(0);
    for (i = 1; i < A.size(); ++i)
    {
        ArrayItem val = A[i];
        if (prev > val) break;
        prev = val;
        A.mark(i);
    }

    if (i == A.size()) {
        return true;
    }

    // unmark
    while (i > 0) A.unmark(i--);
    A.unmark(0);

    return false;
}

// ****************************************************************************

QuickSortPivotType g_quicksort_pivot = PIVOT_FIRST;

// some quicksort variants use hi inclusive and some exclusive, we require it
// to be _exclusive_. hi == array.end()!
ssize_t QuickSortSelectPivot(SortArray& A, ssize_t lo, ssize_t hi)
{
    if (g_quicksort_pivot == PIVOT_FIRST)
        return lo;

    if (g_quicksort_pivot == PIVOT_LAST)
        return hi-1;

    if (g_quicksort_pivot == PIVOT_MID)
        return (lo + hi) / 2;

    if (g_quicksort_pivot == PIVOT_RANDOM)
        return lo + (rand() % (hi - lo));

    if (g_quicksort_pivot == PIVOT_MEDIAN3)
    {
        ssize_t mid = (lo + hi) / 2;

        // cases if two are equal
        if (A[lo] == A[mid]) return lo;
        if (A[lo] == A[hi-1] || A[mid] == A[hi-1]) return hi-1;

        // cases if three are different
        return A[lo] < A[mid]
            ? (A[mid] < A[hi-1] ? mid : (A[lo] < A[hi-1] ? hi-1 : lo))
            : (A[mid] > A[hi-1] ? mid : (A[lo] < A[hi-1] ? lo : hi-1));
    }

    if (g_quicksort_pivot == PIVOT_MOM)
    {
        volatile ssize_t j = hi;
        ssize_t k; 
        A.watch(&j);

        while(j > lo + 1)
        {
            k = j;
            j = lo;
            for(ssize_t i = lo; i < k; i += 5)
            {
                InsertionSortExtra(A, i, std::min(i+5, k));
                A.swap(j++, (i + std::min(i+5, k)) / 2);
            }
        }
        A.unwatch_all();
        
        return lo;
    }

    return lo;
}

wxArrayString QuickSortPivotText()
{
    wxArrayString sl;

    sl.Add( _("First Item") );
    sl.Add( _("Last Item") );
    sl.Add( _("Middle Item") );
    sl.Add( _("Random Item") );
    sl.Add( _("Median of Three") );
    sl.Add( _("Median of Medians") );

    return sl;
}

// ****************************************************************************
// *** Selection Sort

void SelectionSort(SortArray& A, bool stable)
{
    volatile ssize_t jMin = 0;
    A.watch(&jMin, 3);

    for (size_t i = 0; i < A.size()-1; ++i)
    {
        jMin = i;

        for (size_t j = i+1; j < A.size(); ++j)
        {
            if (A[j] < A[jMin]) {
                A.mark_swap(j, jMin);
                jMin = j;
            }
        }

        if(stable) {
            for(int rot = jMin; rot > i; rot--)
            {
                A.swap(rot, rot - 1);
            }
        } else {
            A.swap(i, jMin);
        }

        // mark the last good element
        if (i > 0) A.unmark(i-1);
        A.mark(i);
    }
    A.unwatch_all();
}

void SelectionSort(SortArray& A) { SelectionSort(A, false); }
void StableSelectionSort(SortArray& A) { SelectionSort(A, true); }

// ****************************************************************************
// *** Insertion Sort

// swaps every time (keeps all values visible)
void InsertionSort(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type key = A[i];
        A.mark(i);

        ssize_t j = i - 1;
        while (j >= 0 && A[j] > key)
        {
            A.swap(j, j+1);
            j--;
        }

        A.unmark(i);
    }
}

// swaps every time (keeps all values visible)
void BinaryInsertionSort(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type key = A[i];
        A.mark(i);

        int lo = 0, hi = i;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (key <= A[mid])
                hi = mid;
            else
                lo = mid + 1;
        }

        // item has to go into position lo

        ssize_t j = i - 1;
        while (j >= lo)
        {
            A.swap(j, j+1);
            j--;
        }

        A.unmark(i);
    }
}

// ****************************************************************************
// *** Merge Sort (out-of-place with sentinels)

// by Timo Bingmann

void Merge(SortArray& A, size_t lo, size_t mid, size_t hi)
{
    // mark merge boundaries
    A.mark(lo);
    A.mark(mid,3);
    A.mark(hi-1);

    // allocate output
    std::vector<value_type> out(hi-lo);

    // merge
    size_t i = lo, j = mid, o = 0; // first and second halves
    while (i < mid && j < hi)
    {
        // copy out for fewer time steps
        value_type ai = A[i], aj = A[j];

        out[o++] = (ai < aj ? (++i, ai) : (++j, aj));
    }

    // copy rest
    while (i < mid) out[o++] = A[i++];
    while (j < hi) out[o++] = A[j++];

    ASSERT(o == hi-lo);

    A.unmark(mid);

    // copy back
    for (i = 0; i < hi-lo; ++i)
        A.set(lo + i, out[i]);

    A.unmark(lo);
    A.unmark(hi-1);
}

void MergeSort(SortArray& A, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        size_t mid = (lo + hi) / 2;

        MergeSort(A, lo, mid);
        MergeSort(A, mid, hi);

        Merge(A, lo, mid, hi);
    }
}

void MergeSort(SortArray& A)
{
    return MergeSort(A, 0, A.size());
}

void MergeSortIterative(SortArray& A)
{
    for (size_t s = 1; s < A.size(); s *= 2)
    {
        for (size_t i = 0; i + s < A.size(); i += 2 * s)
        {
            Merge(A, i, i + s,
                  std::min(i + 2 * s, A.size()));
        }
    }
}

// ****************************************************************************
// *** Quick Sort LR (in-place, pointers at left and right, pivot is middle element)

// by Timo Bingmann, based on Hoare's original code

void QuickSortLR(SortArray& A, ssize_t lo, ssize_t hi)
{
    // pick pivot and watch
    volatile ssize_t p = QuickSortSelectPivot(A, lo, hi+1);

    value_type pivot = A[p];
    A.watch(&p, 2);

    volatile ssize_t i = lo, j = hi;
    A.watch(&i, 3);
    A.watch(&j, 3);

    while (i <= j)
    {
        while (A[i] < pivot)
            i++;

        while (A[j] > pivot)
            j--;

        if (i <= j)
        {
            A.swap(i,j);

            // follow pivot if it is swapped
            if (p == i) p = j;
            else if (p == j) p = i;

            i++, j--;
        }
    }

    A.unwatch_all();

    if (lo < j)
        QuickSortLR(A, lo, j);

    if (i < hi)
        QuickSortLR(A, i, hi);
}

void QuickSortLR(SortArray& A)
{
    return QuickSortLR(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Quick Sort LL (in-place, two pointers at left, pivot is first element and moved to right)

// by Timo Bingmann, based on CLRS' 3rd edition

size_t PartitionLL(SortArray& A, size_t lo, size_t hi)
{
    // pick pivot and move to back
    size_t p = QuickSortSelectPivot(A, lo, hi);

    value_type pivot = A[p];
    A.swap(p, hi-1);
    A.mark(hi-1);

    volatile ssize_t i = lo;
    A.watch(&i, 3);

    for (size_t j = lo; j < hi-1; ++j)
    {
        if (A[j] <= pivot) {
            A.swap(i, j);
            ++i;
        }
    }

    A.swap(i, hi-1);
    A.unmark(hi-1);
    A.unwatch_all();

    return i;
}

void QuickSortLL(SortArray& A, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        size_t mid = PartitionLL(A, lo, hi);

        QuickSortLL(A, lo, mid);
        QuickSortLL(A, mid+1, hi);
    }
}

void QuickSortLL(SortArray& A)
{
    return QuickSortLL(A, 0, A.size());
}

// ****************************************************************************
// *** Quick Sort Ternary (in-place, two pointers at left, pivot is first element and moved to right)

// by Timo Bingmann, loosely based on multikey quicksort by B&S

void QuickSortTernaryLR(SortArray& A, ssize_t lo, ssize_t hi)
{
    if (hi <= lo) return;

    int cmp;

    // pick pivot and swap to back
    ssize_t piv = QuickSortSelectPivot(A, lo, hi+1);
    A.swap(piv, hi);
    A.mark(hi);

    const value_type& pivot = A[hi];

    // schema: |p ===  |i <<< | ??? |j >>> |q === |piv
    volatile ssize_t i = lo, j = hi-1;
    volatile ssize_t p = lo, q = hi-1;

    A.watch(&i, 3);
    A.watch(&j, 3);

    for (;;)
    {
        // partition on left
        while (i <= j && (cmp = A[i].cmp(pivot)) <= 0)
        {
            if (cmp == 0) {
                A.mark(p,4);
                A.swap(i, p++);
            }
            ++i;
        }

        // partition on right
        while (i <= j && (cmp = A[j].cmp(pivot)) >= 0)
        {
            if (cmp == 0) {
                A.mark(q,4);
                A.swap(j, q--);
            }
            --j;
        }

        if (i > j) break;

        // swap item between < > regions
        A.swap(i++, j--);
    }

    // swap pivot to right place
    A.swap(i,hi);
    A.mark_swap(i,hi);

    ssize_t num_less = i - p;
    ssize_t num_greater = q - j;

    // swap equal ranges into center, but avoid swapping equal elements
    j = i-1; i = i+1;

    ssize_t pe = lo + std::min(p-lo, num_less);
    for (ssize_t k = lo; k < pe; k++, j--) {
        A.swap(k,j);
        A.mark_swap(k,j);
    }

    ssize_t qe = hi-1 - std::min(hi-1-q, num_greater-1); // one already greater at end
    for (ssize_t k = hi-1; k > qe; k--, i++) {
        A.swap(i,k);
        A.mark_swap(i,k);
    }

    A.unwatch_all();
    A.unmark_all();

    QuickSortTernaryLR(A, lo, lo + num_less - 1);
    QuickSortTernaryLR(A, hi - num_greater + 1, hi);
}

void QuickSortTernaryLR(SortArray& A)
{
    return QuickSortTernaryLR(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Quick Sort LL (in-place, two pointers at left, pivot is first element and moved to right)

// by Timo Bingmann

std::pair<ssize_t,ssize_t> PartitionTernaryLL(SortArray& A, ssize_t lo, ssize_t hi)
{
    // pick pivot and swap to back
    ssize_t p = QuickSortSelectPivot(A, lo, hi);

    value_type pivot = A[p];
    A.swap(p, hi-1);
    A.mark(hi-1);

    volatile ssize_t i = lo, k = hi-1;
    A.watch(&i, 3);

    for (ssize_t j = lo; j < k; ++j)
    {
        int cmp = A[j].cmp(pivot); // ternary comparison
        if (cmp == 0) {
            A.swap(--k, j);
            --j; // reclassify A[j]
            A.mark(k,4);
        }
        else if (cmp < 0) {
            A.swap(i++, j);
        }
    }

    // unwatch i, because the pivot is swapped there
    // in the first step of the following swap loop.
    A.unwatch_all();

    ssize_t j = i + (hi-k);

    for (ssize_t s = 0; s < hi-k; ++s) {
        A.swap(i+s, hi-1-s);
        A.mark_swap(i+s, hi-1-s);
    }
    A.unmark_all();

    return std::make_pair(i,j);
}

void QuickSortTernaryLL(SortArray& A, size_t lo, size_t hi)
{
    if (lo + 1 < hi)
    {
        std::pair<ssize_t,ssize_t> mid = PartitionTernaryLL(A, lo, hi);

        QuickSortTernaryLL(A, lo, mid.first);
        QuickSortTernaryLL(A, mid.second, hi);
    }
}

void QuickSortTernaryLL(SortArray& A)
{
    return QuickSortTernaryLL(A, 0, A.size());
}

// ****************************************************************************
// *** Dual-Pivot Quick Sort

// by Sebastian Wild

void dualPivotYaroslavskiy(class SortArray& a, int left, int right)
{
    if (right > left)
    {
        if (a[left] > a[right]) {
            a.swap(left, right);
        }

        const value_type p = a[left];
        const value_type q = a[right];

        a.mark(left);
        a.mark(right);

        volatile ssize_t l = left + 1;
        volatile ssize_t g = right - 1;
        volatile ssize_t k = l;

        a.watch(&l, 3);
        a.watch(&g, 3);
        a.watch(&k, 3);

        while (k <= g)
        {
            if (a[k] < p) {
                a.swap(k, l);
                ++l;
            }
            else if (a[k] >= q) {
                while (a[g] > q && k < g)  --g;
                a.swap(k, g);
                --g;

                if (a[k] < p) {
                    a.swap(k, l);
                    ++l;
                }
            }
            ++k;
        }
        --l;
        ++g;
        a.swap(left, l);
        a.swap(right, g);

        a.unmark_all();
        a.unwatch_all();

        dualPivotYaroslavskiy(a, left, l - 1);
        dualPivotYaroslavskiy(a, l + 1, g - 1);
        dualPivotYaroslavskiy(a, g + 1, right);
    }
}

void QuickSortDualPivot(class SortArray& a)
{
    return dualPivotYaroslavskiy(a, 0, a.size()-1);
}

// ****************************************************************************
// *** Bubble Sort (Optimized), fixed

void BubbleSort(SortArray& A)
{
    bool swapped;
    for (size_t i = 0; i < A.size()-1; ++i)
    {
        swapped = false;
        for (size_t j = 0; j < A.size()-1 - i; ++j)
        {
            if (A[j] > A[j + 1])
            {
                A.swap(j, j+1);
                swapped = true;
            }
        }
        if(!swapped) { return; }
    }
}

// ****************************************************************************
// *** Cocktail Shaker Sort (More optimized)

// from http://de.wikibooks.org/wiki/Algorithmen_und_Datenstrukturen_in_C/_Shakersort

void CocktailShakerSort(SortArray& A)
{
    size_t lo = 0, hi = A.size()-1, mov = lo;

    while (lo < hi)
    {
        A.mark(lo);
        A.mark(hi);
        for (size_t i = lo; i < hi; ++i)
        {
            if (A[i] > A[i+1])
            {
                A.swap(i, i+1);
                mov = i;
            }
        }
        A.unmark(hi);
        hi = mov;
        A.mark(hi);

        A.mark(lo);
        for (size_t i = hi; i > lo; --i)
        {
            if (A[i-1] > A[i])
            {
                A.swap(i-1, i);
                mov = i;
            }
        }

        A.unmark(lo);
        lo = mov;

    }
}

// ****************************************************************************
// *** Gnome Sort

// from http://en.wikipediA.org/wiki/Gnome_sort

void GnomeSort(SortArray& A)
{
    for (size_t i = 1; i < A.size(); )
    {
        if (A[i] >= A[i-1])
        {
            ++i;
        }
        else
        {
            A.swap(i, i-1);
            if (i > 1) --i;
        }
    }
}

// ****************************************************************************
// *** Comb Sort

// from http://en.wikipediA.org/wiki/Comb_sort

void CombSort(SortArray& A)
{
    const double shrink = 1.3;

    bool swapped = false;
    size_t gap = A.size();

    while ((gap > 1) || swapped)
    {
        if (gap > 1) {
            gap = (size_t)((float)gap / shrink);
        }

        swapped = false;

        for (size_t i = 0; gap + i < A.size(); ++i)
        {
            if (A[i] > A[i + gap])
            {
                A.swap(i, i+gap);
                swapped = true;
            }
        }
    }
}

// ****************************************************************************
// *** Odd-Even Sort

// from http://en.wikipediA.org/wiki/Odd%E2%80%93even_sort

void OddEvenSort(SortArray& A)
{
    bool sorted = false;

    while (!sorted)
    {
        sorted = true;

        for (size_t i = 1; i < A.size()-1; i += 2)
        {
            if(A[i] > A[i+1])
            {
                A.swap(i, i+1);
                sorted = false;
            }
        }

        for (size_t i = 0; i < A.size()-1; i += 2)
        {
            if(A[i] > A[i+1])
            {
                A.swap(i, i+1);
                sorted = false;
            }
        }
    }
}

// ****************************************************************************
// *** Shell Sort

// with gaps by Robert Sedgewick from http://www.cs.princeton.edu/~rs/shell/shell.c

void ShellSort(SortArray& A)
{
    size_t incs[16] = { 1391376, 463792, 198768, 86961, 33936,
                        13776, 4592, 1968, 861, 336,
                        112, 48, 21, 7, 3, 1 };

    for (size_t k = 0; k < 16; k++)
    {
        for (size_t h = incs[k], i = h; i < A.size(); i++)
        {
            value_type v = A[i];
            size_t j = i;

            while (j >= h && A[j-h] > v)
            {
                A.set(j, A[j-h]);
                j -= h;
            }

            A.set(j, v);
        }
    }
}

// ****************************************************************************
// *** Heap Sort

// heavily adapted from http://www.codecodex.com/wiki/Heapsort

bool isPowerOfTwo(size_t x)
{
    return ((x != 0) && !(x & (x - 1)));
}

uint32_t prevPowerOfTwo(uint32_t x)
{
    x |= x >> 1; x |= x >> 2; x |= x >> 4;
    x |= x >> 8; x |= x >> 16;
    return x - (x >> 1);
}

int largestPowerOfTwoLessThan(int n)
{
    int k = 1;
    while (k < n) k = k << 1;
    return k >> 1;
}

void HeapSort(SortArray& A)
{
    size_t n = A.size(), i = n / 2;

    // mark heap levels with different colors
    for (size_t j = i; j < n; ++j)
        A.mark(j, log(prevPowerOfTwo(j+1)) / log(2) + 4);

    while (1)
    {
        if (i > 0) {
            // build heap, sift A[i] down the heap
            i--;
        }
        else {
            // pop largest element from heap: swap front to back, and sift
            // front A[0] down the heap
            n--;
            if (n == 0) return;
            A.swap(0,n);

            A.mark(n);
            if (n+1 < A.size()) A.unmark(n+1);
        }

        size_t parent = i;
        size_t child = i*2 + 1;

        // sift operation - push the value of A[i] down the heap
        while (child < n)
        {
            if (child + 1 < n && A[child + 1] > A[child]) {
                child++;
            }
            if (A[child] > A[parent]) {
                A.swap(parent, child);
                parent = child;
                child = parent*2+1;
            }
            else {
                break;
            }
        }

        // mark heap levels with different colors
        if(n == A.size()) A.mark(i, log(prevPowerOfTwo(i+1)) / log(2) + 4);
    }

}

// ****************************************************************************
// *** Radix Sort (counting sort, most significant digit (MSD) first, in-place redistribute)

// by Timo Bingmann

void RadixSortMSD(SortArray& A, size_t lo, size_t hi, size_t depth)
{
    A.mark(lo); A.mark(hi-1);

    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = floor( log(A.array_max()+1) / log(RADIX) );
    ASSERT(depth <= pmax);

    size_t base = pow(RADIX, pmax - depth);

    // count digits
    std::vector<size_t> count(RADIX, 0);

    for (size_t i = lo; i < hi; ++i)
    {
        size_t r = A[i].get() / base % RADIX;
        ASSERT(r < RADIX);
        count[r]++;
    }

    // inclusive prefix sum
    std::vector<size_t> bkt(RADIX, 0);
    std::partial_sum(count.begin(), count.end(), bkt.begin());

    // mark bucket boundaries
    for (size_t i = 0; i < bkt.size(); ++i) {
        if (bkt[i] == 0) continue;
        A.mark(lo + bkt[i]-1, 3);
    }

    // reorder items in-place by walking cycles
    for (size_t i=0, j; i < (hi-lo); )
    {
        while ( (j = --bkt[ (A[lo+i].get() / base % RADIX) ]) > i )
        {
            A.swap(lo + i, lo + j);
        }
        i += count[ (A[lo+i].get() / base % RADIX) ];
    }

    A.unmark_all();

    // no more depth to sort?
    if (depth+1 > pmax) return;

    // recurse on buckets
    size_t sum = lo;
    for (size_t i = 0; i < RADIX; ++i)
    {
        if (count[i] > 1)
            RadixSortMSD(A, sum, sum+count[i], depth+1);
        sum += count[i];
    }
}

void RadixSortMSD(SortArray& A)
{
    return RadixSortMSD(A, 0, A.size(), 0);
}

// ****************************************************************************
// *** Radix Sort (counting sort, least significant digit (LSD) first, out-of-place redistribute)

// by Timo Bingmann

void RadixSortLSD(SortArray& A)
{
    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = ceil( log(A.array_max()+1) / log(RADIX) );

    for (unsigned int p = 0; p < pmax; ++p)
    {
        size_t base = pow(RADIX, p);

        // count digits and copy data
        std::vector<size_t> count(RADIX, 0);
        std::vector<value_type> copy(A.size());

        for (size_t i = 0; i < A.size(); ++i)
        {
            size_t r = (copy[i] = A[i]).get() / base % RADIX;
            ASSERT(r < RADIX);
            count[r]++;
        }

        // exclusive prefix sum
        std::vector<size_t> bkt(RADIX+1, 0);
        std::partial_sum(count.begin(), count.end(), bkt.begin()+1);

        // mark bucket boundaries
        for (size_t i = 0; i < bkt.size()-1; ++i) {
            if (bkt[i] >= A.size()) continue;
            A.mark(bkt[i], 3);
        }

        // redistribute items back into array (stable)
        for (size_t i=0; i < A.size(); ++i)
        {
            size_t r = copy[i].get() / base % RADIX;
            A.set( bkt[r]++, copy[i] );
        }

        A.unmark_all();
    }
}

// ****************************************************************************
// *** Use STL Sorts via Iterator Adapters

void StlSort(SortArray& A)
{
    std::sort(MyIterator(&A,0), MyIterator(&A,A.size()));
}

void StlStableSort(SortArray& A)
{
    std::stable_sort(MyIterator(&A,0), MyIterator(&A,A.size()));
}

void StlHeapSort(SortArray& A)
{
    std::make_heap(MyIterator(&A,0), MyIterator(&A,A.size()));
    std::sort_heap(MyIterator(&A,0), MyIterator(&A,A.size()));
}

// ****************************************************************************
// *** BogoSort and more slow sorts

// by Timo Bingmann

void BogoSort(SortArray& A)
{
    // keep a permutation of [0,size)
    std::vector<size_t> perm(A.size());

    for (size_t i = 0; i < A.size(); ++i)
        perm[i] = i;

    while (1)
    {
        // check if array is sorted
        if (CheckSorted(A)) break;

        // pick a random permutation of indexes
        std::random_shuffle(perm.begin(), perm.end());

        // permute array in-place
        std::vector<char> pmark(A.size(), 0);

        for (size_t i = 0; i < A.size(); ++i)
        {
            if (pmark[i]) continue;

            // walk a cycle
            size_t j = i;

            //std::cout << "cycle start " << j << " -> " << perm[j] << "\n";

            while ( perm[j] != i )
            {
                ASSERT(!pmark[j]);
                A.swap(j, perm[j]);
                pmark[j] = 1;

                j = perm[j];
                //std::cout << "cycle step " << j << " -> " << perm[j] << "\n";
            }
            //std::cout << "cycle end\n";

            ASSERT(!pmark[j]);
            pmark[j] = 1;
        }

        //std::cout << "permute end\n";

        for (size_t i = 0; i < A.size(); ++i)
            ASSERT(pmark[i]);
    }
}

void BozoSort(SortArray& A)
{
    srand(time(NULL));

    while (1)
    {
        // check if array is sorted
        if (CheckSorted(A)) break;

        // swap two random items
        A.swap(rand() % A.size(), rand() % A.size());
    }
}

// ****************************************************************************
// *** Bitonic Sort as "Parallel" Sorting Network

// from http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm

// modified to first record the recursively generated swap sequence, and then
// sort it back into the order a parallel sorting network would perform the
// swaps in

namespace BitonicSortNetworkNS {

struct swappair_type
{
    // swapped positions
    unsigned int i,j;

    // depth of recursions: sort / merge
    unsigned int sort_depth, merge_depth;

    swappair_type(unsigned int _i, unsigned int _j,
                  unsigned int _sort_depth, unsigned int _merge_depth)
        : i(_i), j(_j),
          sort_depth(_sort_depth), merge_depth(_merge_depth)
    { }

    // order relation for sorting swaps
    bool operator < (const swappair_type& b) const
    {
        if (sort_depth != b.sort_depth)
            return sort_depth > b.sort_depth;

        if (merge_depth != b.merge_depth)
            return merge_depth < b.merge_depth;

        return i < b.i;
    }
};

typedef std::vector<swappair_type> sequence_type;
std::vector<swappair_type> sequence;

void replay(SortArray& A)
{
    for (sequence_type::const_iterator si = sequence.begin();
         si != sequence.end(); ++si)
    {
        if (A[si->i] > A[si->j])
            A.swap(si->i, si->j);
    }
}

static const bool ASCENDING = true; // sorting direction

static void compare(SortArray& /* A */, unsigned int i, unsigned int j, bool dir,
                    unsigned int sort_depth, unsigned int merge_depth)
{
    // if (dir == (A[i] > A[j])) A.swap(i, j);

    if (dir)
        sequence.push_back( swappair_type(i,j, sort_depth, merge_depth) );
    else
        sequence.push_back( swappair_type(j,i, sort_depth, merge_depth) );
}

static void bitonicMerge(SortArray& A, unsigned int lo, unsigned int n, bool dir,
                         unsigned int sort_depth, unsigned int merge_depth)
{
    if (n > 1)
    {
        unsigned int m = largestPowerOfTwoLessThan(n);

        for (unsigned int i = lo; i < lo + n - m; i++)
            compare(A, i, i + m, dir, sort_depth, merge_depth);

        bitonicMerge(A, lo, m, dir, sort_depth, merge_depth+1);
        bitonicMerge(A, lo + m, n - m, dir, sort_depth, merge_depth+1);
    }
}

static void bitonicSort(SortArray& A, unsigned int lo, unsigned int n, bool dir,
                        unsigned int sort_depth)
{
    if (n > 1)
    {
        unsigned int m = n / 2;
        bitonicSort(A, lo, m, !dir, sort_depth+1);
        bitonicSort(A, lo + m, n - m, dir, sort_depth+1);
        bitonicMerge(A, lo, n, dir, sort_depth, 0);
    }
}

void sortIt(SortArray& A)
{
    sequence.clear();
    bitonicSort(A, 0, A.size(), true, 0);
    std::sort(sequence.begin(), sequence.end());
    replay(A);
    sequence.clear();
}

void sortRe(SortArray& A)
{
    sequence.clear();
    bitonicSort(A, 0, A.size(), true, 0);
    replay(A);
    sequence.clear();
}

} // namespace BitonicSortNS

void BitonicSortIt(SortArray& A)
{
    BitonicSortNetworkNS::sortIt(A);
}
void BitonicSortRe(SortArray& A)
{
    BitonicSortNetworkNS::sortRe(A);
}

// ****************************************************************************
// *** Batcher's Odd-Even Merge Sort as "Parallel" Sorting Network

// from http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/networks/oemen.htm

// modified to first record the recursively generated swap sequence, and then
// sort it back into the order a parallel sorting network would perform the
// swaps in

namespace BatcherSortNetworkNS {

struct swappair_type
{
    // swapped positions
    unsigned int i,j;

    // depth of recursions: sort / merge
    unsigned int sort_depth, merge_depth;

    swappair_type(unsigned int _i, unsigned int _j,
                  unsigned int _sort_depth, unsigned int _merge_depth)
        : i(_i), j(_j),
          sort_depth(_sort_depth), merge_depth(_merge_depth)
    { }

    // order relation for sorting swaps
    bool operator < (const swappair_type& b) const
    {
        if (sort_depth != b.sort_depth)
            return sort_depth > b.sort_depth;

        if (merge_depth != b.merge_depth)
            return merge_depth > b.merge_depth;

        return i < b.i;
    }
};

typedef std::vector<swappair_type> sequence_type;
std::vector<swappair_type> sequence;

void replay(SortArray& A)
{
    for (sequence_type::const_iterator si = sequence.begin();
         si != sequence.end(); ++si)
    {
        if (A[si->i] > A[si->j])
            A.swap(si->i, si->j);
    }
}

static void compare(SortArray& A, unsigned int i, unsigned int j,
                    unsigned int sort_depth, unsigned int merge_depth)
{
    // skip all swaps beyond end of array
    ASSERT(i < j);
    if (j >= A.size()) return;

    sequence.push_back( swappair_type(i,j, sort_depth, merge_depth) );

    //if (A[i] > A[j]) A.swap(i, j);
}

// lo is the starting position and n is the length of the piece to be merged, r
// is the distance of the elements to be compared
static void oddEvenMerge(SortArray& A, unsigned int lo, unsigned int n, unsigned int r,
                         unsigned int sort_depth, unsigned int merge_depth)
{
    unsigned int m = r * 2;
    if (m < n)
    {
        // even subsequence
        oddEvenMerge(A, lo, n, m, sort_depth, merge_depth+1);
        // odd subsequence
        oddEvenMerge(A, lo + r, n, m, sort_depth, merge_depth+1);

        for (unsigned int i = lo + r; i + r < lo + n; i += m)
            compare(A, i, i + r, sort_depth, merge_depth);
    }
    else {
        compare(A, lo, lo + r, sort_depth, merge_depth);
    }
}

// sorts a piece of length n of the array starting at position lo
static void oddEvenMergeSort(SortArray& A, unsigned int lo, unsigned int n,
                             unsigned int sort_depth)
{
    if (n > 1)
    {
        unsigned int m = n / 2;
        oddEvenMergeSort(A, lo, m, sort_depth+1);
        oddEvenMergeSort(A, lo + m, m, sort_depth+1);
        oddEvenMerge(A, lo, n, 1, sort_depth, 0);
    }
}

void sortIt(SortArray& A)
{
    sequence.clear();

    unsigned int n = largestPowerOfTwoLessThan(A.size());
    if (n != A.size()) n *= 2;

    oddEvenMergeSort(A, 0, n, 0);
    std::sort(sequence.begin(), sequence.end());
    replay(A);
    sequence.clear();
}

void sortRe(SortArray& A)
{
    sequence.clear();

    unsigned int n = largestPowerOfTwoLessThan(A.size());
    if (n != A.size()) n *= 2;

    oddEvenMergeSort(A, 0, n, 0);
    replay(A);
    sequence.clear();
}

} // namespace BatcherSortNetworkNS

void BatcherSortNetworkIt(SortArray& A)
{
    BatcherSortNetworkNS::sortIt(A);
}

void BatcherSortNetworkRe(SortArray& A)
{
    BatcherSortNetworkNS::sortRe(A);
}

// ****************************************************************************
// *** Smooth Sort

// from http://en.wikipediA.org/wiki/Smoothsort

namespace SmoothSortNS {

static const int LP[] = {
    1, 1, 3, 5, 9, 15, 25, 41, 67, 109,
    177, 287, 465, 753, 1219, 1973, 3193, 5167, 8361, 13529, 21891,
    35421, 57313, 92735, 150049, 242785, 392835, 635621, 1028457,
    1664079, 2692537, 4356617, 7049155, 11405773, 18454929, 29860703,
    48315633, 78176337, 126491971, 204668309, 331160281, 535828591,
    866988873 // the next number is > 31 bits.
};

static void sift(SortArray& A, int pshift, int head)
{
    // we do not use Floyd's improvements to the heapsort sift, because we
    // are not doing what heapsort does - always moving nodes from near
    // the bottom of the tree to the root.

    value_type val = A[head];

    while (pshift > 1)
    {
        int rt = head - 1;
        int lf = head - 1 - LP[pshift - 2];

        if (val.cmp(A[lf]) >= 0 && val.cmp(A[rt]) >= 0)
            break;

        if (A[lf].cmp(A[rt]) >= 0) {
            A.set(head, A[lf]);
            head = lf;
            pshift -= 1;
        }
        else {
            A.set(head, A[rt]);
            head = rt;
            pshift -= 2;
        }
    }

    A.set(head, val); 
}

static void trinkle(SortArray& A, int p, int pshift, int head, bool isTrusty)
{
    value_type val = A[head];

    while (p != 1)
    {
        int stepson = head - LP[pshift];

        if (A[stepson].cmp(val) <= 0)
            break; // current node is greater than head. sift.

        // no need to check this if we know the current node is trusty,
        // because we just checked the head (which is val, in the first
        // iteration)
        if (!isTrusty && pshift > 1) {
            int rt = head - 1;
            int lf = head - 1 - LP[pshift - 2];
            if (A[rt].cmp(A[stepson]) >= 0 ||
                A[lf].cmp(A[stepson]) >= 0)
                break;
        }

        A.set(head, A[stepson]);

        head = stepson;
        //int trail = Integer.numberOfTrailingZeros(p & ~1);
        int trail = __builtin_ctz(p & ~1);
        p >>= trail;
        pshift += trail;
        isTrusty = false;
    }

    if (!isTrusty) {
        A.set(head, val);
        sift(A, pshift, head);
    }
}

void sort(SortArray& A, int lo, int hi)
{
    int head = lo; // the offset of the first element of the prefix into m

    // These variables need a little explaining. If our string of heaps
    // is of length 38, then the heaps will be of size 25+9+3+1, which are
    // Leonardo numbers 6, 4, 2, 1.
    // Turning this into a binary number, we get b01010110 = 0x56. We represent
    // this number as a pair of numbers by right-shifting all the zeros and
    // storing the mantissa and exponent as "p" and "pshift".
    // This is handy, because the exponent is the index into L[] giving the
    // size of the rightmost heap, and because we can instantly find out if
    // the rightmost two heaps are consecutive Leonardo numbers by checking
    // (p&3)==3

    int p = 1; // the bitmap of the current standard concatenation >> pshift
    int pshift = 1;

    while (head < hi)
    {
        if ((p & 3) == 3) {
            // Add 1 by merging the first two blocks into a larger one.
            // The next Leonardo number is one bigger.
            sift(A, pshift, head);
            p >>= 2;
            pshift += 2;
        }
        else {
            // adding a new block of length 1
            if (LP[pshift - 1] >= hi - head) {
                // this block is its final size.
                trinkle(A, p, pshift, head, false);
            } else {
                // this block will get merged. Just make it trusty.
                sift(A, pshift, head);
            }

            if (pshift == 1) {
                // LP[1] is being used, so we add use LP[0]
                p <<= 1;
                pshift--;
            } else {
                // shift out to position 1, add LP[1]
                p <<= (pshift - 1);
                pshift = 1;
            }
        }
        p |= 1;
        head++;
    }

    trinkle(A, p, pshift, head, false);

    while (pshift != 1 || p != 1)
    {
        if (pshift <= 1) {
            // block of length 1. No fiddling needed
            //int trail = Integer.numberOfTrailingZeros(p & ~1);
            int trail = __builtin_ctz(p & ~1);
            p >>= trail;
            pshift += trail;
        }
        else {
            p <<= 2;
            p ^= 7;
            pshift -= 2;

            // This block gets broken into three bits. The rightmost bit is a
            // block of length 1. The left hand part is split into two, a block
            // of length LP[pshift+1] and one of LP[pshift].  Both these two
            // are appropriately heapified, but the root nodes are not
            // necessarily in order. We therefore semitrinkle both of them

            trinkle(A, p >> 1, pshift + 1, head - LP[pshift] - 1, true);
            trinkle(A, p, pshift, head - 1, true);
        }

        head--;
    }
}

} // namespace SmoothSortNS

void SmoothSort(SortArray& A)
{
    return SmoothSortNS::sort(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Stooge Sort

void StoogeSort(SortArray& A, int i, int j, bool stable)
{
    if ((j - i + 1 <= 2 || !stable) && A[i] > A[j])
    {
        A.swap(i, j);
    }

    if (j - i + 1 >= 3)
    {
        int t = (j - i + 1) / 3;

        A.mark(i, 3);
        A.mark(j, 3);

        StoogeSort(A, i, j-t, stable);
        StoogeSort(A, i+t, j, stable);
        StoogeSort(A, i, j-t, stable);

        A.unmark(i);
        A.unmark(j);
    }
}

void StoogeSort(SortArray& A)
{
    StoogeSort(A, 0, A.size()-1, false);
}

void StableStoogeSort(SortArray& A)
{
    StoogeSort(A, 0, A.size()-1, true);
}

// ****************************************************************************
// *** Slow Sort

void SlowSort(SortArray& A, int i, int j)
{
    if (i >= j) return;

    int m = (i + j) / 2;

    SlowSort(A, i, m);
    SlowSort(A, m+1, j);

    if (A[m] > A[j])
        A.swap(m, j);

    A.mark(j, 2);

    SlowSort(A, i, j-1);

    A.unmark(j);
}

void SlowSort(SortArray& A)
{
    SlowSort(A, 0, A.size()-1);
}

// ****************************************************************************
// *** Cycle Sort

// Adapted from http://en.wikipedia.org/wiki/Cycle_sort

void CycleSort(SortArray& array, ssize_t n)
{
    volatile ssize_t cycleStart = 0;
    array.watch(&cycleStart, 16);

    volatile ssize_t rank = 0;
    array.watch(&rank, 3);

    // Loop through the array to find cycles to rotate.
    for (cycleStart = 0; cycleStart < n - 1; ++cycleStart)
    {
        value_type& item = array.get_mutable(cycleStart);

        do {
            // Find where to put the item.
            rank = cycleStart;
            for (ssize_t i = cycleStart + 1; i < n; ++i)
            {
                if (array[i] < item)
                    rank++;
            }

            // If the item is already there, this is a 1-cycle.
            if (rank == cycleStart) {
                array.mark(rank, 2);
                break;
            }

            // Otherwise, put the item after any duplicates.
            while (item == array[rank])
                rank++;

            // Put item into right place and colorize
            std::swap(array.get_mutable(rank), item);
            array.mark(rank, 2);

            // Continue for rest of the cycle.
        }
        while (rank != cycleStart);
    }

    array.unwatch_all();
}

void CycleSort(SortArray& A)
{
    CycleSort(A, A.size());
}


// ******************************************************************************
// ** Some More Sorts
// ******************************************************************************

// ****************************************************************************
// ** Circle and Circloid

// Circloid by yuji

bool CircleSorting(SortArray& A, int lo, int hi, bool mode) 
{ 
    A.mark(hi,3);

    bool swapped = false;
    if (hi > lo) {
        int m = (hi - lo + 1)/2;

        // Circloid version
        if(mode) {
            swapped = CircleSorting(A, lo, hi-m, true) || swapped;
            swapped = CircleSorting(A, lo+m, hi, true) || swapped;
        }

        for(int off = 0;off < m; ++off) {
            if (A[lo+off] > A[hi-off]) { A.swap(lo+off, hi-off); swapped = true; }
        }

        // Circle version
        if(!mode) {
            swapped = CircleSorting(A, lo, hi-m, false) || swapped;
            swapped = CircleSorting(A, lo+m, hi, false) || swapped;
        }
    }

    A.unmark(hi);

    return swapped;
} 

void CircleSort(SortArray& A)
{
    while( CircleSorting(A, 0, A.size() - 1, false) );
}

void CircloidSort(SortArray& A) // custom sort.
{
    while( CircleSorting(A, 0, A.size() - 1, true) ); 
}

// ********************************************************************************
// ** Pairwise Sorting Network as "Parallel" Sorting Network + Recursive version

// implemented by Piotr, translation to C++ by me (u-ndefined)

// modified to first record the recursively generated swap sequence, and then
// sort it back into the order a parallel sorting network would perform the
// swaps in

namespace PairwiseSortingNetworkNS {

struct swappair_type
{
    // swapped positions
    unsigned int i,j;

    // depth of recursions: sort / merge
    unsigned int sort_depth, merge_depth;

    swappair_type(unsigned int _i, unsigned int _j,
                  unsigned int _sort_depth, unsigned int _merge_depth)
        : i(_i), j(_j),
          sort_depth(_sort_depth), merge_depth(_merge_depth)
    { }

    // order relation for sorting swaps
    bool operator < (const swappair_type& b) const
    {
        if (sort_depth != b.sort_depth)
            return sort_depth > b.sort_depth;

        if (merge_depth != b.merge_depth)
            return merge_depth > b.merge_depth;

        return i < b.i;
    }
    

};

typedef std::vector<swappair_type> sequence_type;
std::vector<swappair_type> sequence;

void replay(SortArray& A)
{
    for (sequence_type::const_iterator si = sequence.begin();
         si != sequence.end(); ++si)
    {
        if (A[si->i] > A[si->j])
            A.swap(si->i, si->j);
    }
}

static void compare(SortArray& A, unsigned int i, unsigned int j,
                    unsigned int sort_depth, unsigned int merge_depth)
{
    // skip all swaps beyond end of array
    ASSERT(i < j);
    if (j >= A.size()) return;

    sequence.push_back( swappair_type(i,j, sort_depth, merge_depth) );

    //if (A[i] > A[j]) A.swap(i, j);
}

void pairwiserecursive(SortArray& A, int start, int end, int gap,
    unsigned int sort_depth) {

     if (start == end - gap){
        return;
    }

    int b = start + gap;

     while (b < end){
        compare(A, b - gap, b, -sort_depth-1, A.size());
        b += (2 * gap);
    }
    if (((end - start) / gap) % 2 == 0) {
        pairwiserecursive(A, start, end, gap * 2, sort_depth+1);
        pairwiserecursive(A, start + gap, end + gap, gap * 2, sort_depth+1);
    } else {
        pairwiserecursive(A, start, end + gap, gap * 2, sort_depth+1);
        pairwiserecursive(A, start + gap, end, gap * 2, sort_depth+1);
    }

    int a = 1;
    while (a < ((end - start) / gap)) {
        a = (a * 2) + 1;
    }
    b = start + gap;
    while (b + gap < end){
        int c = a;
        while (c > 1) {
            c /= 2;
            if (b + (c * gap) < end){
                compare(A, b, b + (c * gap), sort_depth, c * gap);
            }
        }
        b += (2 * gap);
    }
}

void sortIt(SortArray& A)
{
    sequence.clear();
    pairwiserecursive(A, 0, A.size(), 1, 0);
    std::sort(sequence.begin(), sequence.end());
    replay(A);
    sequence.clear();
}

void sortRe(SortArray& A)
{
    sequence.clear();
    pairwiserecursive(A, 0, A.size(), 1, 0);
    replay(A);
    sequence.clear();
}

} // namespace PairwiseSortingNetworkNS

void PairwiseSortingIt(SortArray& A) {
    PairwiseSortingNetworkNS::sortIt(A);
}

void PairwiseSortingRe(SortArray& A) {
    PairwiseSortingNetworkNS::sortRe(A);
}

// ****************************************************************************
// ** Stable Quick Sort

void StableQuickSort(SortArray& A, size_t lo, size_t hi)
{
    if((hi - lo) > 1) {
        size_t j = QuickSortSelectPivot(A,lo,hi);
        volatile ssize_t loa_end = lo;
        A.watch(&loa_end, 3);

        std::vector<value_type> loa;
        std::vector<value_type> hia;
        value_type pivot = A[j];
        for(size_t i = lo; i < hi; i++)
        {
            if(A[i] < pivot) {
               loa.push_back(A[i]);
               loa_end++;
            } else {
               hia.push_back(A[i]);
            }
        }

        int p = loa_end;
        j = 0;
        while(j < loa.size())
        {
            A.set(j + lo, loa[j]);
            j++;
        }

        A.unwatch_all();

        while((j - p + lo) < hia.size())
        {
            A.set(j + lo, hia[j - p + lo]);
            j++;
        }
        hia.clear();
        loa.clear();

        StableQuickSort(A, lo, p);
        StableQuickSort(A, p + 1, hi);
    }
}

void StableQuickSort(SortArray& A)
{
    StableQuickSort(A, 0, A.size());
}

// ********************************************************************************
// ** Bose-Nelson Sorting Network as "Parallel" Sorting Network + Recursive version

// implemented by atinm, modified by me (u-ndefined)

// modified to first record the recursively generated swap sequence, and then
// sort it back into the order a parallel sorting network would perform the
// swaps in

namespace BoseNelsonNetworkNS {

struct swappair_type
{
    // swapped positions
    unsigned int i,j;

    // depth of recursions: sort / merge
    unsigned int sort_depth, merge_depth;

    swappair_type(unsigned int _i, unsigned int _j,
                  unsigned int _sort_depth, unsigned int _merge_depth)
        : i(_i), j(_j),
          sort_depth(_sort_depth), merge_depth(_merge_depth)
    { }

    // order relation for sorting swaps
    bool operator < (const swappair_type& b) const
    {
        if (sort_depth != b.sort_depth)
            return sort_depth > b.sort_depth;

        if (merge_depth != b.merge_depth)
            return merge_depth > b.merge_depth;

        return i < b.i;
    }
};

typedef std::vector<swappair_type> sequence_type;
std::vector<swappair_type> sequence;

void replay(SortArray& A)
{
    for (sequence_type::const_iterator si = sequence.begin();
         si != sequence.end(); ++si)
    {
        if (A[si->i] > A[si->j])
            A.swap(si->i, si->j);
    }
}

static void compare(SortArray& A, unsigned int i, unsigned int j,
                    unsigned int sort_depth, unsigned int merge_depth)
{
    // skip all swaps beyond end of array
    ASSERT(i < j);
    if (j >= A.size()) return;

    sequence.push_back( swappair_type(i,j, sort_depth, merge_depth) );

    //if (A[i] > A[j]) A.swap(i, j);
}

void bosemerge(SortArray& A, unsigned int i,  /* value of first element in sequence 1 */
         unsigned int x,  /* length of sequence 1 */
         unsigned int j,  /* value of first element in sequence 2 */
         unsigned int y,  /* length of sequence 2 */
         unsigned int sort_depth)  
{
    unsigned int a, b;

    if(x == 1 && y == 1) { compare(A, i, j, sort_depth, j - i); } /* 1 comparison sorts 2 items */
    else if(x == 1 && y == 2)
    {
        /* 2 comparisons inserts an item into an
         * already sorted sequence of length 2. */
        compare(A, i, j + 1, sort_depth, (j + 1) - i);
        compare(A, i, j, sort_depth, j - i);
    }
    else if(x == 2 && y == 1)
    {
        compare(A, i, j, sort_depth, j - i);
        compare(A, i + 1, j, sort_depth, j - (i + 1));
    }
    else
    {
        /* Recurse on shorter sequences, attempting
         * to make the length of one subsequence odd
         * and the length of the other even. If we
         * can do this, we eventually merge the two. */
        a = x/2;
        b = (x & 1) ? (y/2) : ((y + 1)/2);
        bosemerge(A, i, a, j, b, sort_depth);
        bosemerge(A, (i + a), (x - a), (j + b), (y - b), sort_depth);
        bosemerge(A, (i + a), (x - a), j, b, sort_depth);
    }
}

void bosepartition(SortArray& A, unsigned int i, /* value of first element in sequence */ unsigned int m, /* length of sequence */
    unsigned int sort_depth) 
{
    int a;

    if(m > 1)
    {
        /* Partition into 2 shorter sequences,
         * generate a sorting method for each,
         * and merge the two sub-networks. */
        a = m/2;
        bosepartition(A, i, a, sort_depth+1);
        bosepartition(A, (i + a), (m - a), sort_depth+1);
        bosemerge(A, i, a, (i + a), (m - a), sort_depth);
    }
}

void sortIt(SortArray& A)
{
    sequence.clear();
    bosepartition(A, 0, A.size(), 0);
    std::sort(sequence.begin(), sequence.end());
    replay(A);
    sequence.clear();
}

void sortRe(SortArray& A)
{
    sequence.clear();
    bosepartition(A, 0, A.size(), 0);
    replay(A);
    sequence.clear();
}


} // namespace BoseNelsonNetworkNS

void BoseNelsonSortingIt(SortArray& A)
{
    BoseNelsonNetworkNS::sortIt(A); /* sort the sequence {X1,...,Xn} */
}

void BoseNelsonSortingRe(SortArray& A)
{
    BoseNelsonNetworkNS::sortRe(A); /* sort the sequence {X1,...,Xn} */
}

// ****************************************************************************
// ** Gravity Sort

void GravitySort(SortArray& A)
{
    int sz = A.size();
    int max = A.array_max();

    // Make abacus
    bool* abacus = new bool[sz*max];
    // Fill abacus
	for (int i = 0; i < sz; i++) {
		int tmp = A[i].get();
		for (int j = 0; j < tmp; j++)
			abacus[i*max+j] = true;
		// Fill rest with zeroes
		for (int j = tmp; j < max; j++) {
			abacus[i*max+j] = false;
		}
	}

    // simulate abacus
	for (int i = max - 1; i >= 0; i--) {
		volatile ssize_t sum = 0;
		for (int j = 0; j < sz; j++) {
			if (abacus[i+j*max]) {
				sum++;
                int it = A[j].get();
                A.set(j,ArrayItem(it-1));
			}
		}
        A.mark(sz - sum);
		for (int j = sz - 1; j >= sz - sum; j--) {
            int it = A[j].get();
            A.set(j,ArrayItem(it+1));
		}
        A.unmark(sz - sum);
	}
}

// ****************************************************************************
// ** Less Bogo Sort

void LessBogoSort(SortArray& A)
{
    volatile ssize_t cur = 0;
    A.watch(&cur);

    // keep a permutation of [0,size)
    std::vector<size_t> perm(A.size());
    for (size_t i = 0; i < A.size(); ++i)
        perm[i] = i;

    while (cur < A.size() - 1)
    {
        size_t k = cur;
        while(++k < A.size()) {
            if(A[k] < A[cur]) {break;}
        }
        if(k == A.size()) {
            cur++;
            perm.pop_back();
            for (size_t i = 0; i < A.size() - cur; ++i)
                perm[i] = i;
        } else {
            // pick a random permutation of indexes
            std::random_shuffle(perm.begin(), perm.end());
            // permute array in-place
            std::vector<char> pmark(A.size() - cur, 0);

            for (size_t i = 0; i < A.size() - cur; ++i)
            {
                if (pmark[i]) continue;

                // walk cycle

                size_t j = i;
                while ( perm[j] != i )
                {
                    ASSERT(!pmark[j]);
                    A.swap(j + cur, perm[j] + cur);
                    pmark[j] = 1;

                    j = perm[j];
                }

                ASSERT(!pmark[j]);
                pmark[j] = 1;
            }

            for (size_t i = 0; i < A.size() - cur; ++i)
                ASSERT(pmark[i]);
        }
    }

    A.unwatch_all();
}

// ****************************************************************************
// ** Flipped Min Heap Sort

void FlippedMinHeapSort(SortArray& A)
{
    size_t n = A.size(), i = n / 2;
    ssize_t k = 0;

    // mark heap levels with different colors
    for (size_t j = 0; j < i; ++j) A.mark(j, log(prevPowerOfTwo(n-j)) / log(2) + 4);

    i--;

    while (1)
    {
        if (i < (n - 1)) {
            // build heap, sift A[i] down the heap
            i++;
        }
        else {
        
            // pop smallest element from heap: swap front to back, and sift
            // front A[0] down the heap
            if (k == n) return;
            A.swap(k,n-1);

            A.mark(k);
            if ((k-1) >= 0) A.unmark(k-1);
            k++;
        }

        ssize_t parent = i;
        ssize_t child = n - ((n-i)*2);


        // sift operation - push the value of A[i] down the heap
        while (child >= k)
        {
            if (((child - 1) >= k) && (A[child - 1] < A[child])) {
                child--;
            }
            if (A[child] < A[parent]) {
                A.swap(parent, child);
                parent = child;
                child = n - ((n - parent)*2);
            }
            else {
                break;
            }
        }


        // mark heap levels with different colors
        A.mark(i, log(prevPowerOfTwo(n-i)) / log(2) + 4);
    }

}

// ****************************************************************************
// *** Iterative Circle Sort

// uses a queue to keep track of the ranges

bool CircleSortingIt(SortArray& A, int lo, int hi)
{
    bool swapped = false;
    std::queue<int> ranges;

    // push ranges
    ranges.push(lo);
    ranges.push(hi);
    while(!ranges.empty())
    {
        // get ranges
        int qlo = ranges.front();
        ranges.pop();
        int qhi = ranges.front();
        ranges.pop();

        int m = (qhi - qlo + 1)/2;
        for (int off = 0; off < m; ++off)
        {
            if (A[qlo+off] > A[qhi-off]) { A.swap(qlo+off, qhi-off); swapped = true; }
        }

        if(qhi - qlo >= 2)
        {
            ranges.push(qlo);
            ranges.push(qhi-m);

            ranges.push(qlo+m);
            ranges.push(qhi);
        }
    }
    return swapped;
}

void IterativeCircleSort(SortArray& A) {
    while( CircleSortingIt(A, 0, A.size( ) - 1) );
}

// ****************************************************************************
// *** Balanced Sorting Network

void BalancedSortingNetwork(SortArray& A) {
    int k = 0;
    for (int num = ceil(log(A.size() - 1) / log(2.0)); k < num; ++k) {
        CircleSortingIt(A, 0, A.size( ) - 1);
    }
}

// ****************************************************************************
// *** Proxmap Sort

void ProxmapSort(SortArray& A)
{
    std::vector<size_t> hits(A.size());
    std::vector<size_t> maps(A.size());
    std::vector<int> A2(A.size(), -1);

    // compute the array's minimum and maximum
    value_type min = A[0], max = A[0];

    for(int i = 1; i < A.size(); ++i)
    {
        if(A[i] < min)
        {
            min = A[i];
        } else if(A[i] > max)
        {
            max = A[i];
        }
    }

    int imin = min.get(), imax = max.get();

    // compute hits and key maps
    for(int i = 0; i < A.size(); ++i)
    {
        maps[i] = int(float(A[i].get() - imin) / float(imax - imin) * float(A.size() - 1));
        ++hits[maps[i]];
    }

    // prefix sum - using the vector hits[] to compute the proxmaps
    hits[A.size() - 1] = A.size() - hits[A.size() - 1];

    for(int i = A.size() - 1; i > 0; --i)
    {
        hits[i - 1] = hits[i] - hits[i - 1];
    }

    // insert A[i] to A2 in correct position
    int iIdx, iLo;
    for(int i = 0; i < A.size(); ++i)
    {
        iIdx = hits[maps[i]];
        iLo = iIdx;
        int cur = A[i].get();
        while(A2[iIdx] != -1) iIdx++;
        while(iIdx > iLo && cur < A2[iIdx - 1])
        {
            A.touch(iIdx, 1, 1);
            A2[iIdx] = A2[iIdx - 1];
            iIdx--;
        }
        A.touch(iIdx, 1, 1);
        A2[iIdx] = cur;
    }
    for(int i = 0; i < A.size(); ++i) A.set(i, (ArrayItem)A2[i]); 
}

/*/ ****************************************************************************
// *** Improved Weave Merge Sort 

// in-shuffle based on https://arxiv.org/pdf/0805.1598.pdf

void ImprovedWeaveMergeSort(SortArray& A, size_t lo, size_t hi)
{
    int l = hi - lo;
    if(l < 2) return;
    size_t m = (lo + hi) / 2;
    // ImprovedWeaveMergeSort(A, lo, m);
    // ImprovedWeaveMergeSort(A, m, hi);

    // the in-shuffle algorithm

    int st = 0;
    while(st < l)
    {
        // step 1: find m
        int i=3;
        while(i <= l+1) i*=3;

        if(i > 3) i=i/3;

        // step 2: rotation (NOTE: currently uses reversals)

        // step 3: do the cycle leader algorithm 
        int idx, tmp1, tmp2;

        idx = (st*2)%(l+1);
        tmp1 = A[idx-1];
        A[idx-1] = A[st-1];

        while(idx != st)
        {
            tmp2 = A[(idx * 2)%(l + 1)-1];
            A[(idx * 2)%(l + 1)-1] = tmp1;
            tmp1 = tmp2;
            idx = (idx * 2) % (l+1);
        }

        // step 4: recursive repeat
    }

    InsertionSortExtra(A, lo, hi);
}

void ImprovedWeaveMergeSort(SortArray& A)
{
    ImprovedWeaveMergeSort(A, 0, A.size());
}

// ****************************************************************************/

namespace UnbalancedTreeSortNS {

// BST node that keeps track of array indices for visualization
struct Node
{
    size_t parent;
    size_t depth;
    Node* left;
    Node* right;
};

// special types of function for visualization
// indexInsert is like insert function of BST but for indices
Node* indexInsert(SortArray& A, Node* parent, size_t i, size_t depth)
{
    if(parent == NULL)
    {
        parent = new Node();
        parent->parent = i;
        parent->depth = depth;
        parent->left = NULL;
        parent->right = NULL;
    } else if(A[i] >= A[parent->parent]) {
        parent->right = indexInsert(A, parent->right, i, depth+1);
    } else {
        parent->left = indexInsert(A, parent->left, i, depth+1);
    }

    return parent;
}

// likewise, this is an in-order traversal function of BST but modified to work for indices
// traversal - left subtree - parent - right subtree
size_t traverse(SortArray& A, std::vector<ArrayItem> aux, Node* parent, size_t i)
{

    if(parent->left != NULL)
    {
        i = traverse(A, aux, parent->left, i);
    }
    A.set(i++, aux[parent->parent]);
    if(parent->right != NULL)
    {
        i = traverse(A, aux, parent->right, i);
    }

    // used to retrieve the number of counts for the array
    return i;
}


void sort(SortArray& A)
{
    Node* tree;
    tree = NULL;
    std::vector<ArrayItem> aux(A.size());

    for(size_t i=0; i<A.size(); i++)
    {
        tree = indexInsert(A, tree, i, 0);
        aux[i] = A[i]; // copy to another array
    }
    traverse(A, aux, tree, 0);
}

} // namespace UnbalancedTreeSortNS

void UnbalancedTreeSort(SortArray& A)
{
    UnbalancedTreeSortNS::sort(A);
}