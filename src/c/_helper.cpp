//
// Created by Colin Small on 7/1/21.
//

#include <vector>
#include <string>
#include <cmath>



/**
 * Helper function to return multiple items from a string given a list of indices to access.
 * @param str The string to access.
 * @param indices A vector containing the indices to access and return from str.
 * @return A new vector containing only the elements at the specified indices.
 */
std::vector<char> accessMultipleIndicesString(std::string str, std::vector<int>& indices)
{
    // Make sure indices isn't larger than the array
    assert(str.size() >= indices.size());
    // Make sure max index in indices isn't larger than length of array
    assert(*std::max_element(indices.begin(), indices.end()) < str.size());

    std::vector<char> returnVec = {};
    returnVec.reserve(indices.size());

    for(int i = 0; i < indices.size(); i++){
        returnVec.push_back(str[i]);
    }

    return returnVec;
}

/**
 * Returns a vector of indices where characters in a string are equal to a given character
 * @param str The string to get indices from.
 * @param comparator The character to be compared to the contents of str.
 * @return A vector of indices.
 */
std::vector<int> getIndicesWhereEqualFromString(std::string_view str, char comparator){
    std::vector<int> returnVec = {};

    for(int i = 0; i < str.size(); i++){
        if(str[i] == comparator)
            returnVec.push_back(i);
    }

    return returnVec;
}