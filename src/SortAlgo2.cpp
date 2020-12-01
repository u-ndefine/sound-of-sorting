/******************************************************************************
 * src/SortAlgo2.cpp
 *
 * Most of the custom (user-made) sorting algorithms are implemented here.
 * https://www.youtube.com/playlist?list=PL5w_-zMAJC8tGIdBbUCoBFkAkiVMQClhZ 
 *
 * Some algorithms are annotated using the mark() and watch()
 * functions from SortArray.
 *
 ******************************************************************************
 * The algorithms in this file are copyrighted by the original authors. All
 * code is freely available.
 *
 * All modifications are copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
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

// ********************************************************************************
// ** Y - Slow Sort

// by yuji

void YSlowSort(SortArray& arr, int lo, int hi) 
{ 
    if (arr[lo] > arr[hi]) { arr.swap(lo, hi); }
    if (hi - lo > 1) {
        int m = (hi - lo + 1)/2.0;
        for (int run = 0; run <= 1; run++) 
        {
            
            arr.mark(hi, 3);

            YSlowSort(arr, lo, hi-m);
            YSlowSort(arr, lo+m, hi);

            arr.unmark(hi);

            arr.mark(lo, 2);
            arr.mark(hi, 2);

            YSlowSort(arr, lo+1, hi-1);
            
            arr.unmark(lo);
            arr.unmark(hi);
        }
    }
} 

void YSlowSort(SortArray& A)
{
    YSlowSort(A, 0, A.size() - 1);
}

// ****************************************************************************
// ** Quad Stooge Sort

// sort idea by EilrahcF

void QuadStoogeSort(SortArray& A, size_t lo, size_t len)
{
    if (len >= 2 && A[lo] > A[lo+len-1])
    {
        A.swap(lo, lo+len-1);
    }

    if (len >= 3)
    {
		size_t len1 = len/2;
		size_t len2 = (len+1)/2;
		size_t len3 = (len1+1)/2 + (len2+1)/2;
		
        A.mark(lo+len1, 3);
		QuadStoogeSort(A, lo, len1);
        A.unmark(lo+len1);
		QuadStoogeSort(A, lo+len1, len2);
		QuadStoogeSort(A, lo+len1/2, len3);
        A.mark(lo+len1, 3);
		QuadStoogeSort(A, lo+len1, len2);
        A.unmark(lo+len1);
		QuadStoogeSort(A, lo, len1);
		if(len > 3) //at len == 3 this makes a useless comparison
			QuadStoogeSort(A, lo+len1/2, len3);
    }
}

void QuadStoogeSort(SortArray& A)
{
    QuadStoogeSort(A, 0, A.size());
}

// ****************************************************************************
// ** Bubble Scan Quicksort

// by aphitorite + thatsOven

void BubbleScanQuicksort(SortArray& A, int lo, int hi)
{
    while((hi - lo) > 16) {

        // bubble select method
        int sum = 0;
        bool swapped = false;
        for(int i = lo; i < hi - 1; i++) {
            if (A[i] > A[i+1]){ A.swap(i, i+1); swapped = true; }
            sum += A[i].get();
        }

        if(!swapped){ return; }
        // partition
        int mean = (int)(sum / float(hi - lo - 1));

        volatile ssize_t i = lo, j = hi - 1;
        A.watch(&i, 3);
        A.watch(&j, 3);

        while (i <= j)
        {
            while (A[i].get() <= mean && i <= j)
                i++;

            while (A[j].get() > mean && i <= j)
                j--;

            if (i <= j) { A.swap(i++,j--); }
        }

        A.unwatch_all();
        BubbleScanQuicksort(A, (int)(j + 1), hi - 1);
        hi = j + 1;
    }

    for (size_t i = lo; i < hi; ++i)
    {
        value_type key = A[i];

        ssize_t j = i - 1;
        while (j >= lo && A[j] > key)
        {
            A.swap(j, j+1);
            j--;
        }
    }
}

void BubbleScanQuicksort(SortArray& A)
{
    BubbleScanQuicksort(A, 0, A.size());
}

// ****************************************************************************
// ** Kiyomi, Kinnow, Chinotto, and Mandarin Orange Sort

// by yuji

void GappedCocktailSort(SortArray& A, bool mode1, bool mode2)
{
    bool swapped = true;
    int gap = 1;
    while(swapped)
    {
        swapped = false;
        int i = 0;
        while((i + gap) < A.size())
        {
            if(A[i] > A[i + gap])
            {
                A.mark(i + gap);
                A.mark(i);

                if(mode1) {
                    for(int rot = i; rot < i + gap; rot++)
                    {
                        A.swap(rot, rot + 1);
                    }
                } else {
                    for(int rot = i + gap; rot > i; rot--)
                    {
                        A.swap(rot, rot - 1);
                    }
                }
                A.unmark(i + gap);
                A.unmark(i);

                swapped = true;
                gap += 1;
            } else if (gap >= 2) {
                gap -= 1;
            }
            i++;
        }
        while((i - gap) > 0)
        {
            if(A[i] < A[i - gap])
            {
                A.mark(i - gap);
                A.mark(i);

                if(mode2) {
                    for(int rot = i - gap; rot < i; rot++)
                    {
                        A.swap(rot, rot + 1);
                    }
                } else {
                    for(int rot = i; rot > i - gap; rot--)
                    {
                        A.swap(rot, rot - 1);
                    }
                }
                A.unmark(i - gap);
                A.unmark(i);
                swapped = true;
                gap += 1;
            } else if (gap >= 2) {
                gap -= 1;
            }
            i--;
        }
    }
}

void KiyomiSort(SortArray& A)
{
    GappedCocktailSort(A, false, true);

    // failsafe
    while(!CheckSorted(A)) {
        GappedCocktailSort(A, false, true);
    }
}
void KinnowSort(SortArray& A)
{
    GappedCocktailSort(A, true, true);
}
void ChinottoSort(SortArray& A)
{
    GappedCocktailSort(A, true, false);
}
void MandarinOrangeSort(SortArray& A)
{
    GappedCocktailSort(A, false, false);
    
    // failsafe
    while(!CheckSorted(A)) {
        GappedCocktailSort(A, false, false);
    }
}

// ****************************************************************************
// ** Optimized Stable Quick Sort

// by aphitorite

void OptiStableQuickSort(SortArray& A, size_t lo, size_t hi)
{
    if((hi - lo) > 1) {

        size_t j = QuickSortSelectPivot(A,lo,hi);
        volatile ssize_t loa_end = lo;
        A.watch(&loa_end, 3);

        std::vector<value_type> aux;
        int pivot = A[j].get();
        for(size_t i = lo; i < hi; i++)
        {
            ++g_compare_count; // count as comparison
            if(A[i].get() < pivot) {
               A.set(loa_end++, A[i]);
            } else {
               aux.push_back(A[i]);
            }
        }

        int p = loa_end;

        A.unwatch_all();

        while((loa_end - p) < aux.size())
        {
            A.set(loa_end, aux[loa_end - p]);
            loa_end++;
        }
        aux.clear();

        OptiStableQuickSort(A, lo, p);
        OptiStableQuickSort(A, p + 1, hi);
    }
}

void OptiStableQuickSort(SortArray& A)
{
    OptiStableQuickSort(A, 0, A.size());
}

// ****************************************************************************
// ** Float Sort

// by Lance, implemented by me (u-ndefined)

void FloatSort(SortArray& A)
{
    bool swapped = true;
    while(swapped)
    {
        swapped = false;
        for(int cur = 0; cur < A.size() - 1; cur++)
        {
            A.mark(cur);
            int i = cur - 1;
            while(i > 0 && A[i-1] > A[i]) {
                A.swap(i-1,i);
                i--;
                swapped = true;
            }
            while(i < A.size() - 1 && A[i] > A[i+1]) {
                A.swap(i,i+1);
                i++;
                swapped = true;
            }
            A.unmark(cur);
        }
    }
}

// ****************************************************************************
// ** Quicksort (LLL pointers)

// by myself (u-ndefined)

void QuickSortLLL(SortArray& A, size_t lo, size_t hi)
{
    if ((hi - lo) > 1) {
        // pick pivots and move them to back
        size_t p1 = QuickSortSelectPivot(A, lo, hi);
        A.swap(p1, hi-1);
    
        size_t p2 = QuickSortSelectPivot(A, lo, hi - 1);
        A.swap(p2, hi-2);

        if(A[hi-2] > A[hi-1]) { A.swap(hi-1,hi-2); }

        A.mark(hi-1);
        A.mark(hi-2);

        volatile ssize_t i = lo;
        A.watch(&i, 3);
        volatile ssize_t j = lo;
        A.watch(&j, 3);

        for (size_t k = lo; k < hi - 1; ++k)
        {
            if (A[k] <= A[hi-1]) {
                A.swap(j, k);
                if(A[j] <= A[hi-2]) {
                    A.swap(i, j);
                   ++i;
                }
                ++j;
            }
        }

        A.swap(hi-1, j);
        A.unmark_all();
        A.unwatch_all();

        QuickSortLLL(A, lo, i-1);
        QuickSortLLL(A, i, j);
        QuickSortLLL(A, j+1, hi);
    }
}

void QuickSortLLL(SortArray& A)
{
    QuickSortLLL(A, 0, A.size());
}

// ****************************************************************************
// ** Grate Sorts

// sort idea by EilrahcF

bool GrateSorting(SortArray& A, int dir) { // returns if swapped
    bool swapped = false;
    volatile ssize_t i;
    if(dir > 0)
    {
        A.watch(&i, 3);
    } else {
        A.watch(&i);
    }
    for(i = 0; i < A.size(); i++)
    {
        size_t j;
        if(dir > 0)
        {
            j = i;
        } else {
            j = A.size();
        }
        j += dir;
        while(j < A.size() && j > i && A[i] <= A[j])
            j += dir;
        if(j < A.size() && j > i)
        {
            A.swap(i, j);
            swapped = true;
        }
    }
    A.unwatch_all();
    return swapped;
}

void GrateSort(SortArray& A) {
    while(GrateSorting(A, -1)) {}
}

void ReverseGrateSort(SortArray& A) {
    while(GrateSorting(A, 1)) {}
}

void CocktailGrateSort(SortArray& A) {
    while(GrateSorting(A, -1)) {
        if(!GrateSorting(A, 1)) {break;}
    }
}

// ****************************************************************************
// ** Jump Down Sort

// by fungamer2

void JumpDownSort(SortArray& A)
{
    volatile ssize_t jMin = A.size() - 1;
    A.watch(&jMin, 3);

    for (size_t i = A.size()-1; i > 0; --i)
    {
        jMin = i;

        for (size_t j = 0; j < i; ++j)
        {
            if (A[j] > A[jMin]) {
                A.swap(j, jMin);
            }
        }

        // mark the last good element
        if (i < A.size()-1) A.unmark(i+1);
        A.mark(i);
    }
    A.unwatch_all();
}

// ****************************************************************************
// ** Room Sort

// by EilrahcF

void RoomSort(SortArray& A)
{
    size_t sz = (size_t)sqrt(A.size());
    volatile ssize_t cur = 1;
    volatile ssize_t limit = 0;
    A.watch(&cur);
    A.watch(&limit);

    for (ssize_t left = A.size(); left > 0; left -= sz)
    {
        A.mark(left-1, 3);
        for (size_t i = 1; i < left; ++i)
        {
            cur = i;
            limit = std::max((ssize_t)0, cur - (ssize_t)sz);
            value_type tmp, key = A[i];

            ssize_t j = i - 1;
            while (j >= limit && (tmp = A[j]) > key)
            {
                A.set(j + 1, tmp);
                j--;
            }
            A.set(j + 1, key);

        }
        A.unmark(left-1);
    }
    A.unwatch_all();
}

// ****************************************************************************
// *** Slope Sort

// by EilrahcF

void SlopeSort(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type key = A[i];
        A.mark(i);

        for(size_t j = i; j > 0; --j)
        {
            if(key < A[j - 1]) A.swap(j, j-1);
        }

        A.unmark(i);
    }
}


// ****************************************************************************
// ** Pancake Sorts

void PancakeSort(SortArray& A)
{
    volatile ssize_t jMax = 0;
    A.watch(&jMax, 3);
    
    for (size_t i = A.size(); i > 0; --i)
    {
        jMax = 0;

        for(size_t j = 0; j < i; ++j)
        {
            if(A[jMax] < A[j])
            {
                jMax = j;
            }
        }

        if(A[i-1] != A[jMax]) {
            std::reverse(MyIterator(&A,0), MyIterator(&A,(size_t)jMax+1));

            if(i < A.size()) A.unmark(i);
            A.mark(i-1);
        
            std::reverse(MyIterator(&A,0), MyIterator(&A,i));
        }

    }

    A.unwatch_all();
}

void QuasiPancakeSort(SortArray& A)
{
    for (size_t i = 1; i < A.size(); ++i)
    {
        value_type key = A[i];
        A.mark(i);

        ssize_t j = i - 1;
        while (j >= 0 && A[j] > key)
        {
            j--;
        }


        std::reverse(MyIterator(&A,j+1), MyIterator(&A,i));
        std::reverse(MyIterator(&A,j+1), MyIterator(&A,i+1));

        A.unmark(i);
    }
}

void BinaryQuasiPancakeSort(SortArray& A)
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

        std::reverse(MyIterator(&A,hi), MyIterator(&A,i));
        std::reverse(MyIterator(&A,hi), MyIterator(&A,i+1));
        
        A.unmark(i);

    }
}

// ****************************************************************************
// *** Bingo Sort

// a variant of selection sort for equal items

void BingoSort(SortArray& A)
{
    size_t i = A.size();
    value_type curMax;
    bool chk = false;

    volatile ssize_t jMaxP = 0;
    A.watch(&jMaxP, 3);

    while (i > 0)
    {
        jMaxP = 0;
        A.mark(i-1);

        for (size_t j = 0; j < i; ++j)
        {
            if (chk && A[j] == curMax) {
                A.swap(--i, j--);
                if(i > 0) {
                    A.mark(i-1);
                }
                A.mark(i,4);
            } else if (A[j] > A[jMaxP]) {
                jMaxP = j;
            }
        }

        chk = true;
        curMax = A[jMaxP];
        A.unmark_all();
    }
    A.unwatch_all();
}

// ****************************************************************************
// *** Iterative Stooge Sort

// uses a queue to retrieve ranges

void IterativeStoogeSort(SortArray& A)
{
    std::queue<int> ranges;

    ranges.push(0);
    ranges.push(A.size() - 1);

    while(!ranges.empty()) // while queue is not empty
    {
        // get ranges
        int i = ranges.front();
        ranges.pop();
        int j = ranges.front();
        ranges.pop();

        // compare and swap
        if(A[i] > A[j]) {A.swap(i, j);}

        if(j - i >= 2)
        {
            int m = (j - i) / 3;
            
            // push the next ranges
            ranges.push(i);
            ranges.push(j-m);

            ranges.push(i+m);
            ranges.push(j);

            ranges.push(i);
            ranges.push(j-m);
        }
    }
}

// ****************************************************************************
// *** Hyperstooge Sort

// by fungamer2

void HyperStoogeSort(SortArray& A, int i, int j)
{
    if (A[i] > A[j])
    {
        A.swap(i, j);
    }

    if (j - i >= 2)
    {

        HyperStoogeSort(A, i, j-1);
        HyperStoogeSort(A, i+1, j);

        A.mark(j);
        HyperStoogeSort(A, i, j-1);
        A.unmark(j);
    }
}

void HyperStoogeSort(SortArray& A)
{
    HyperStoogeSort(A, 0, A.size() - 1);
}

// ****************************************************************************
// *** Reverse Insertion Sort

// with extra item on stack

void ReverseInsertionSort(SortArray& A)
{
    for (ssize_t i = A.size() - 2; i >= 0; --i)
    {
        value_type tmp, key = A[i];
        A.mark(i);

        ssize_t j = i + 1;
        while (j < A.size() && (tmp = A[j]) < key)
        {
            A.set(j - 1, tmp);
            j++;
        }
        A.set(j - 1, key);

        A.unmark(i);
    }
}

// ****************************************************************************
// *** Balance Sort

// sort idea by EilrahcF

void BalanceSort(SortArray& A)
{
    for(size_t i = A.size(); i > 0; --i)
    {
        A.mark(i - 1, 2);

        size_t mid = (i + 1) / 2;
        volatile ssize_t pntL = 0;
        volatile ssize_t pntR = mid;

        A.watch(&pntL, 3);
        A.watch(&pntR, 3);

        A.mark(mid, 2);

        int k;
        while(pntL <= mid && pntR < i)
        {
            if(A[pntL] >= A[pntR]) {
                k = pntL;
                pntR++;
            } else {
                k = pntR;
                pntL++;
            }
        }
        A.swap(i - 1, k);
        A.unmark_all();
        A.unwatch_all();
    }

    InsertionSortExtra(A, 0, A.size());
}

// ****************************************************************************
// *** Triangular Heap Sort

// sort idea by EilrahcF

void TriangularHeapSort(SortArray& A)
{
    size_t n = A.size(), i = 0, szT = 0;

    // find largest triangular number smaller than n
    while(i + szT < n)
    {
        ++szT;
        i += szT;
    }
    // use szT as number of elements not sifted down
    // initialize previous largest triangular number smaller than n
    size_t smT = i - szT;
    i = n - szT;

    for(size_t j = n - 1; j >= i; --j) { //mark different colors of triangular heap
        if((smT + szT) >= j) {
            A.mark(j, 4+(szT % 12));
        } else {
            A.mark(j, 4+((szT + 1) % 12));
        }
    }

    while (1)
    {
        if (i > 0) {
            // build heap, sift A[i] down the heap
            i--;
            if (i < smT) { // check if number lowers the triangular number
                --szT;
                smT -= szT;
            }
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


        size_t k = szT;
        size_t parent = i;
        size_t child = i + k;

        // sift operation - push the value of A[i] down the heap
        while (child < n)
        {
            if (child + 1 < n && A[child + 1] > A[child]) {
                child++;
            }
            if (A[child] > A[parent]) {
                A.swap(parent, child);
                parent = child;
                k++;
                child = parent+k;
            }
            else {
                break;
            }
        }

        // mark different colors of triangular heap
        A.mark(i, 4+(szT % 12));
    }
}

// ****************************************************************************
// *** Slide Sort

// by EilrahcF

void SlideSort(SortArray& A)
{
    size_t fdiv = A.size() / 4,
    mid = A.size() / 2,
    tdiv = A.size() * 3 / 4;

    // make 4 sorted lists
    InsertionSortExtra(A, 0, fdiv);
    InsertionSortExtra(A, fdiv, mid);
    InsertionSortExtra(A, mid, tdiv);
    InsertionSortExtra(A, tdiv, A.size());
    
    // make subarrays
    std::vector<value_type> sub1(fdiv);
    std::vector<value_type> sub2(mid - fdiv);
    std::vector<value_type> sub3(tdiv - mid);
    std::vector<value_type> sub4(A.size() - tdiv);

    for(size_t i = 0; i < fdiv; ++i) sub1[i] = A[i];
    for(size_t i = fdiv; i < mid; ++i) sub2[i - fdiv] = A[i];
    for(size_t i = mid; i < tdiv; ++i) sub3[i - mid] = A[i];
    for(size_t i = tdiv; i < A.size(); ++i) sub4[i - tdiv] = A[i];

    // quad merge
    // this can probably be more optimized
    for (int a = A.size() - 1; a >= 0; --a) {
    value_type largest;
    if (sub1.size() >= 1) {
        largest = sub1.back();
        if (sub2.size() >= 1 && largest < sub2.back()) {
            largest = sub2.back();
        }
        if (sub3.size() >= 1 && largest < sub3.back()) {
            largest = sub3.back();
        }
        if (sub4.size() >= 1 && largest < sub4.back()) {
            largest = sub4.back();
        }
    }
    else if (sub2.size() >= 1) {
        largest = sub2.back();
        if (sub3.size() >= 1 && largest < sub3.back()) {
            largest = sub3.back();
        }
        if (sub4.size() >= 1 && largest < sub4.back()) {
            largest = sub4.back();
        }
    }
    else if (sub3.size() >= 1) {
        largest = sub3.back();
        if (sub4.size() >= 1 && largest < sub4.back()) {
            largest = sub4.back();
        }
    }
    else {
        largest = sub4.back();
    }
    A.set(a, largest);
    if (sub1.size() >= 1 && largest == sub1.back()) {
        sub1.pop_back();
    }
    else if (sub2.size() >= 1 && largest == sub2.back()) {
        sub2.pop_back();
    }
    else if (sub3.size() >= 1 && largest == sub3.back()) {
        sub3.pop_back();
    }
    else {
        sub4.pop_back();
    }
    }
}

// ****************************************************************************
// *** Iterative Quick Sort (LL pointers)

void IterativeQuickSortLL(SortArray& A)
{
    std::vector<ssize_t> stack(A.size()*2);
    int top = -1;

    // push range
    stack[++top] = 0;
    stack[++top] = A.size();

    while(top >= 0) // while stack is not empty
    {
        // get ranges
        ssize_t hi = stack[top--];
        ssize_t lo = stack[top--];

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

        // push partitions
        if(i > lo) {
            stack[++top] = lo;
            stack[++top] = i;
        }
        if(i + 1 < hi) {
            stack[++top] = i + 1;
            stack[++top] = hi;
        }
    }
}

// ****************************************************************************
// *** Iterative Quick Sort (LR pointers)

void IterativeQuickSortLR(SortArray& A)
{
    std::vector<ssize_t> stack(A.size()*2);
    int top = -1;

    // push range
    stack[++top] = 0;
    stack[++top] = A.size() - 1;

    while(top >= 0) // while stack is not empty
    {
        // get ranges
        ssize_t hi = stack[top--];
        ssize_t lo = stack[top--];
        
        // partition
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

        // push partitions
        if(lo < j) {
            stack[++top] = lo;
            stack[++top] = j;
        }
        if(i < hi) {
            stack[++top] = i;
            stack[++top] = hi;
        }
    }
}

// ****************************************************************************
// *** Shove Sort

// by Piotr

void ShoveSort(SortArray& A)
{
    // check if sorted
    size_t i, j;
    while(1) {
    value_type prev = A[0];
    A.mark(0);
    for (i = 1; i < A.size(); ++i)
    {
        value_type val = A[i];
        if (prev > val) break;
        prev = val;
        A.mark(i);
    }

    if (i == A.size()) {
        return;
    }

    // unmark and move element to last
    j = i - 1;
    while (i > 0) A.unmark(i--);
    A.unmark(0);

    for(j; j+1 < A.size(); ++j) A.swap(j, j+1);
    }
}

// ****************************************************************************
// *** Iterative Linked Quick Sort (LL pointers)

// uses a queue to retrieve the ranges

void IterativeLinkedQuickSortLL(SortArray& A)
{
    std::queue<ssize_t> ranges;

    // queue range
    ranges.push(0);
    ranges.push(A.size());

    while(!ranges.empty()) // while queue is not empty
    {
        // get ranges
        ssize_t lo = ranges.front();
        ranges.pop();
        ssize_t hi = ranges.front();
        ranges.pop();

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

        // push partitions
        if(i > lo) {
            ranges.push(lo);
            ranges.push((size_t)i);
        }
        if(i + 1 < hi) {
            ranges.push((size_t)i + 1);
            ranges.push(hi);
        }
    }
}

// ****************************************************************************
// *** Iterative Linked Quick Sort (LR pointers)

// uses a queue to retrieve the ranges

void IterativeLinkedQuickSortLR(SortArray& A)
{
    std::queue<ssize_t> ranges;

    // queue range
    ranges.push(0);
    ranges.push(A.size() - 1);

    while(!ranges.empty()) // while queue is not empty
    {
        // get ranges
        ssize_t lo = ranges.front();
        ranges.pop();
        ssize_t hi = ranges.front();
        ranges.pop();

        // partition
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

        // push partitions
        if(lo < j) {
            ranges.push(lo);
            ranges.push((size_t)j);
        }
        if(i < hi) {
            ranges.push((size_t)i);
            ranges.push(hi);
        }
    }
}

// ****************************************************************************
// *** Rotate Radix Sorts

// by aphitorite

namespace RotateRadixSortNS {

void rotate(SortArray& A, int lo, int sz, int hi) {
    while (sz < hi && sz > lo) {
        if (hi - sz < sz - lo) {
            for (int i = sz; i < hi; ++i) {
                A.swap(i, i - (hi - sz));
            }
            int tmp = hi;
            hi = sz;
            sz -= tmp - sz;
        }
        else {
            for (int i = lo; i < sz; ++i) {
                A.swap(i, i + (sz - lo));
            }
            int tmp = lo;
            lo = sz;
            sz += sz - tmp;
        }
    }
}

void rotateRec(SortArray& A,
    int posA, int posB,
    std::vector<int> a, std::vector<int> b,
    int lo, int hi)
{
    int mid = (lo + hi) / 2;
    if (mid == lo) return;

    rotate(A, posA + a[mid], posB + b[lo], posB + b[mid]);

    int lenA = posB + b[lo] - (posA + a[mid]);
    int lenB = b[mid] - b[lo];
    rotateRec(A, posA + lenB, posB, a, b, mid, hi);
    rotateRec(A, posA, posB - lenA, a, b, lo, mid);
}

std::vector<int> rotateRecMain(SortArray& A,
    int a, int b,
    unsigned int base, const unsigned int RADIX)
{
    if (b - a <= 1) {
        int i = A[a].get() / base % RADIX + 1;
        A.touch(a, 16, 1, 10);
        std::vector<int> regs(RADIX);
        while (i < RADIX) {
            regs[i++] = 1;
        }
        return regs;
    }
    int m = (a + b) / 2;
    std::vector<int> c = rotateRecMain(A, a, m, base, RADIX);
    std::vector<int> d = rotateRecMain(A, m, b, base, RADIX);
    rotateRec(A, a, m, c, d, 0, RADIX);
    for (int j = 0; j < c.size(); ++j) {
        c[j] += d[j];
    }
    return c;
}

void rotateMSD(SortArray& A, int lo, int hi, int p, const unsigned int RADIX) {
    if (hi - lo < 1 || p < 0) {
        return;
    }
    std::vector<int> pos = rotateRecMain(A, lo, hi, pow(RADIX, p), RADIX);
    for (int i = 1; i < pos.size(); ++i) {
        rotateMSD(A, lo + pos[i - 1], lo + pos[i], p - 1, RADIX);
    }
    rotateMSD(A, lo + pos.back(), hi, p - 1, RADIX);
}

} // namespace RotateRadixSortNS

void RotateRadixLSD(SortArray& A) {
    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = floor( log(A.array_max()+1) / log(RADIX) );

    for (unsigned int p = 0; p <= pmax; ++p)
    {
        size_t base = pow(RADIX, p);
        RotateRadixSortNS::rotateRecMain(A, 0, A.size(), base, RADIX);
    }
}

void RotateRadixMSD(SortArray& A) {
    // radix and base calculations
    const unsigned int RADIX = 4;

    unsigned int pmax = floor( log(A.array_max()+1) / log(RADIX) );

    RotateRadixSortNS::rotateMSD(A, 0, A.size(), pmax, RADIX);
}

// ****************************************************************************
// *** Wiggle Sort

// by EilrahcF

void WiggleSort(SortArray& A, size_t lo, size_t hi)
{
    if(hi - lo <= 1) return;
    size_t m = (lo + hi) / 2;
    volatile ssize_t j = m;
    int dir = 1;
    A.watch(&j, 3);
    A.mark(m);
    A.mark(hi-1);
    for(size_t i = lo; i < m; i++)
    {
        while((dir == 1 && j < hi) || (dir == -1 && j >= m))
        {
            if(A[i] > A[j]) A.swap(i,j);
            j += dir;
        }
        dir = -dir;
        j += dir;
    }
    A.unmark_all();
    A.unwatch_all();

    WiggleSort(A, lo, m);
    WiggleSort(A, m, hi);
}

void WiggleSort(SortArray& A)
{
    WiggleSort(A, 0, A.size());
}

// ****************************************************************************
// *** Bogo Bogo Bogo Sort

// by EilrahcF

void BogoBogoBogoSort(SortArray& A)
{
    // keep a permutation of [0,size)
    std::vector<size_t> perm(A.size());
    std::vector<value_type> match(A.size());

    for (size_t i = 0; i < A.size(); ++i)
    {
        match[i] = A[i];
        perm[i] = i;
    }

    while (1)
    {
        std::random_shuffle(match.begin(), match.end());

        while(1) {

            size_t pChk;
            // check if array matches with the permuation
            for(pChk = 0; pChk < A.size(); ++pChk)
            {
                A.mark(pChk, 16);
                if(A[pChk] != match[pChk])
                {
                    break;
                }
            }
            A.unmark_all();
            if(pChk == A.size()) break;

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

        // check if array is sorted
        if (CheckSorted(A)) break;
    }
}
 
// ****************************************************************************
// *** Odd Even Sort (Base 3)

// by Lance

void OddEvenBase3(SortArray& A)
{
    int off = 0;
    int chk = 3;

    while (--chk > 0)
    {
        for (size_t i = off; i < A.size()-1; i += 3)
        {
            if(A[i] > A[i+1])
            {
                A.swap(i, i+1);
                chk = 3;
            }
        }
        off = (off + 1) % 3;
    }
}

// ****************************************************************************
// *** Odd Even Bogo Sort 

// by fungamer2

void OddEvenBogoSort(SortArray& A)
{
    int off = 0;

    while (!CheckSorted(A))
    {
        for (size_t i = off; i < A.size()-1; i += 2)
        {
            size_t ran = (int)(rand() % ((A.size() - off) / 2) ) * 2 + off;
            if(A[ran] > A[ran+1])
            {
                A.swap(ran, ran+1);
            }
        }
        off = 1 - off;
    }
}

// ****************************************************************************
// *** Awkward Sort

// by aphitorite

void AwkwardSort(SortArray& A, int lo, int hi)
{
    if (hi - lo > 1)
    {
        int m = (lo + hi) / 2, l = hi - lo;
        A.mark(hi-1, 3);

        AwkwardSort(A, lo, m);
        AwkwardSort(A, m, hi);
        // compare the halves
        for(int i = lo; (i + l/2) < hi; ++i)
        {
           if(A[i] > A[i + l/2]) {A.swap(i, i+l/2);}
        }
        A.unmark(hi-1);

        AwkwardSort(A, lo + l/4,  lo + l * 3/4); // using an efficient LERP formula
        
        A.mark(hi-1, 3);
        AwkwardSort(A, lo, m);
        AwkwardSort(A, m, hi);
        A.unmark(hi-1);
    }
}

void AwkwardSort(SortArray& A)
{
    AwkwardSort(A, 0, A.size());
}

// ****************************************************************************
// *** Checkerboard Heap Sort

// by EilrahcF

void CheckerboardHeapSort(SortArray& A)
{
    for(size_t i = 1; i < A.size(); ++i)
    {
        size_t parent = i, child = 0;
        while(parent > 0)
        {
            if(parent % 3 == 0)
            {
                child = parent - 2;
                if(child + 1 < A.size() && A[child] < A[child + 1]) child += 1;

            } else {
                child = parent / 3 * 3;
            }
            if(A[parent] < A[child]) A.swap(parent, child);
            else break;
            parent = child;
        }
    }

    for(size_t i = 1; i < A.size() - 1; i += 3)
    {
        if(A[i] > A[i+1]) A.swap(i, i+1);
    }
}

// ****************************************************************************
// *** Modulo Sort

// by McDude_73

void ModuloSort(SortArray& A)
{
    const unsigned int BASE = 2;
    size_t sz = A.size();
    size_t modulo = BASE;
    bool chk = BASE != sz;
    
    std::vector<value_type> aux(A.size());

    while(modulo <= sz + 1)
    {
        for (size_t i = 0; i < sz; ++i) {
            aux[i] = A[i];
        }
        int i = 0;
        int j = 0;
        int k = -1;
        while (j < modulo) {
            int it = aux[i].get();
            if(k+1 < A.size()) A.set(k+1, (ArrayItem)it);
            if (it % modulo == j) ++k;
            if (i + 1 >= sz) {
                ++j;
                i = -1;
            }
            ++i;
        }
        if (modulo * 2 >= sz && chk) {
            chk = false;
            modulo = sz + 1;
        }
        else {
            modulo *= 2;
        }
    }
}

// ****************************************************************************
// *** Cocktail Shell Sort

// by fungamer2

void CocktailShellSort(SortArray& A)
{
    int gap = A.size() / 2;
    bool dir = true; // true -> right, false -> left
    while(gap > 0)
    {
        if(dir) {
            for (ssize_t h = gap, i = h; i < A.size(); i++)
            {
                value_type v = A[i];
                ssize_t j = i;

                while (j >= h && A[j-h] > v)
                {
                    A.set(j, A[j-h]);
                    j -= h;
                }

                A.set(j, v);
            }
        } else {
            for (ssize_t h = gap, i = A.size() - h; i >= 0; i--)
            {
                value_type v = A[i];
                ssize_t j = i;

                while (j < A.size() - h && A[j+h] < v)
                {
                    A.set(j, A[j+h]);
                    j += h;
                }

                A.set(j, v);
            }
        }
        gap = gap / sqrt(2);
        dir = !dir || gap == 1;
    }
}

// ****************************************************************************
// *** Odd-Even Transposition Merge

// by aphitorite

void OddEvenTransMergeSort(SortArray& A, int lo, int hi)
{
    if(hi - lo < 2) return;

    int m = (lo + hi) / 2;
    bool sorted = false;
    OddEvenTransMergeSort(A, lo, m);
    OddEvenTransMergeSort(A, m, hi);

    A.mark(lo);
    A.mark(hi-1);
    for(int i = m - 1; i >= lo && !sorted; --i)
    {
        sorted = true;
        for(int j = i; j < hi - (i - lo) - 1; j += 2)
        {
            if(A[j] > A[j+1])
            {
                A.swap(j, j+1);
                sorted = false;
            }
        }
    }
    for(int i = lo; i <= m && !sorted; ++i)
    {
        sorted = true;
        for(int j = i + 1; j < hi - (i - lo) - 1; j += 2)
        {
            if(A[j] > A[j+1])
            {
                A.swap(j, j+1);
                sorted = false;
            }
        }
    }
    A.unmark(lo);
    A.unmark(hi-1);
}

void OddEvenTransMergeSort(SortArray& A)
{
    OddEvenTransMergeSort(A, 0, A.size());
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

// ****************************************************************************