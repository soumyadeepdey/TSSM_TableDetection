


#include "sorting.h"



	 
/*--------------------------------------------------------Sorting-----------------------------------------------------------------------------------*/

template <typename T>
void IITkgp_functions::swaping(std::vector<T> & data, int i, int j)
{
    T tmp = data[i];
    data[i] = data[j];
    data[j] = tmp;
}


/*-------------------------------------------------------Merge Sort----------------------------------------------------------*/


template <typename T>
void IITkgp_functions::Merge(std::vector<T> & data, int lowl, int highl, int lowr, int highr)
{
    int tmp_low = lowl;
    std::vector<T> tmp;
    
    while (lowl <= highl && lowr <= highr)
    {
        if (data[lowl] < data[lowr])
        {
            tmp.push_back(data[lowl++]);
        }
        else if (data[lowr] < data[lowl])
        {
            tmp.push_back(data[lowr++]);
        }
        else
        {
            tmp.push_back(data[lowl++]);
            tmp.push_back(data[lowr++]);
        }
    }

    while (lowl <= highl)
    {
        tmp.push_back(data[lowl++]);
    }

    while (lowr <= highr)
    {
        tmp.push_back(data[lowr++]);
    }

    
    vector<int>::iterator iter;
    for(iter = tmp.begin(); iter != tmp.end(); ++iter)
    {
        data[tmp_low++] = *iter;
    }
}


template <typename T>
void IITkgp_functions::MergeSort(std::vector<T> & data, int low, int high)
{
    if (low >= high)
    {
        return;
    }
    
    int mid = low + (high-low)/2;

    MergeSort(data, low, mid);

    MergeSort(data, mid+1, high);

    Merge(data, low, mid, mid+1, high);
}



/*-------------------------------------------------------Quick Sort----------------------------------------------------------*/

template <typename T>
int IITkgp_functions::Partition(std::vector<T> & data, int low, int high)
{
    int p = low;
    for (int i = p+1; i <= high; ++i)
    {
        if (data[i] < data[p])
        {
            swaping(data, i, p);
            if (i != p+1)
            {
                swaping(data, i, p+1);
            }
            p = p + 1;
        }
    }

    return p;
}

template <typename T>
void IITkgp_functions::QuickSort(std::vector<T> & data, int low, int high)
{
    if (low >= high) return;
    
    int p = Partition(data, low, high);

    QuickSort(data, low, p-1);
    QuickSort(data, p+1, high);
}




/*-------------------------------------------------------Bubble Sort----------------------------------------------------------*/



//useful for small lists, and for large lists where data is
//already sorted

template <typename T>
void IITkgp_functions::BubbleSort(std::vector<T> & data)
{
    int length = data.size();

    for (int i = 0; i < length; ++i)
    {
        bool swapped = false;
        for (int j = 0; j < length - (i+1); ++j)
        {
            if (data[j] > data[j+1])
            {
                swaping(data, j, j+1);
                swapped = true;
            }
        }
        
        if (!swapped) break;
    }
}

/*-------------------------------------------------------Selection Sort----------------------------------------------------------*/



//useful for small lists and where swapping is expensive
// does at most n swaps
template <typename T>
void IITkgp_functions::SelectionSort(std::vector<T> & data)
{
    int length = data.size();

    for (int i = 0; i < length; ++i)
    {
        int min = i;
        for (int j = i+1; j < length; ++j)
        {
            if (data[j] < data[min])
            {
                min = j;
            }
        }

        if (min != i)
        {
            swaping(data, i, min);
        }
    }
}


/*-------------------------------------------------------Insertion Sort----------------------------------------------------------*/



//useful for small and mostly sorted lists
//expensive to move array elements
template <typename T>
void IITkgp_functions::InsertionSort(std::vector<T> & data)
{
    int length = data.size();

    for (int i = 1; i < length; ++i)
    {
        bool inplace = true;
        int j = 0;
        for (; j < i; ++j)
        {
            if (data[i] < data[j])
            {
                inplace = false;
                break;
            }
        }

        if (!inplace)
        {
            T save = data[i];
            for (int k = i; k > j; --k)
            {
                data[k] = data[k-1];
            }
            data[j] = save;
        }
    }
}

