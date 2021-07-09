//
// Created by Colin Small on 7/1/21.
//

#ifndef CFSNV__HELPER_H
#define CFSNV__HELPER_H

#include <cmath>

template <class C>
int sumMap(std::map<C, int> m){
    int sum = 0;
    typename std::map<C, int>::iterator it;
    for(it = m.begin(); it != m.end(); it++){
        sum += it->second;
    }
    return sum;
}

/**
 * Helper function to return multiple items from a vector given a list of indices to access.
 * @param vec The vector to access.
 * @param indices A vector containing the indices to access and return from array.
 * @return A new vector containing only the elements at the specified indices.
 */
template <class T>
std::vector<T> accessMultipleIndicesVector(std::vector<T>& vec, std::vector<int> indices)
{
//    // Make sure indices isn't larger than the array
//    assert(vec.size() >= indices.size());
//    // Make sure max index in indices isn't larger than length of array
//    assert(*std::max_element(indices.begin(), indices.end()) < vec.size());

    std::vector<T> returnArray(indices.size());

    for(int i = 0; i < indices.size(); i++){
        returnArray[i] = vec[i];
    }

    return returnArray;
}

/**
 * Returns a vector of indices where elements in a vector are equal to some other object
 * TODO: Make this take a generic boolean function so that it can also return indices where items are <, >, etc.
 * TODO: Make this take any type of iterable function.
 * @tparam T The type of objects in the vector and the type that is being compared to
 * @param vec The vector being searched through.
 * @param comparator The object which objects in vec are being compared against.
 * @return
 */
template <class T>
std::vector<int> getIndicesWhereEqualFromVector(std::vector<T>& vec, T comparator){

    std::vector<int> returnVec = {};

    for(int i = 0; i < vec.size(); i++){
        if (vec[i] == comparator)
            returnVec.push_back(i);
    }

    return returnVec;
}

template <class T>
std::vector<int> getIndicesWhereGreaterThan(std::vector<T>& vec, T comparator){

    std::vector<int> returnVec = {};

    for(int i = 0; i < vec.size(); i++){
        if (vec[i] > comparator)
            returnVec.push_back(i);
    }

    return returnVec;
}

// From https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
//    assert(a.size() == b.size());

    std::vector<T> result(a.size(), 0);

//    std::transform(a.begin(), a.end(), b.begin(),
//                   std::back_inserter(result), std::plus<T>());

    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] + b[i];
    }
    return result;
}


template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
//    assert(a.size() == b.size());

    std::vector<T> result(a.size(), 0);

//    std::transform(a.begin(), a.end(), b.begin(),
//                   std::back_inserter(result), std::minus<T>());

    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] - b[i];
    }

    return result;
}

template <typename T>
std::vector<T> operator*(const std::vector<T>& a, const std::vector<T>& b)
{
//    assert(a.size() == b.size());

    std::vector<T> result(a.size(), 0);

//    std::transform(a.begin(), a.end(), b.begin(),
//                   std::back_inserter(result), std::multiplies<T>());

    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] * b[i];
    }
    return result;
}

template <typename T>
std::vector<T> veclog(const std::vector<T>& a)
{
    std::vector<double> loglist(a.size(), 0);

    for(int i = 0; i < a.size(); i ++){
        loglist[i] = (std::log(a[i]));
    }
    return loglist;
}

template <class T>
bool vecContains(std::vector<T> vec, T t){
    return std::find(vec.begin(), vec.end(), t) != vec.end();
}


double sumVector(const std::vector<double>& vec);

double meanVector(const std::vector<double>& vec);

int factorial(int n);

std::vector<std::string> split(const std::string& s, const std::string& delimiter);

std::string toUpper(std::string s);

std::vector<int> getIndicesWhereEqualFromString(std::string_view str, char comparator);
std::vector<char> accessMultipleIndicesString(std::string str, std::vector<int>& indices);


#endif //CFSNV__HELPER_H
