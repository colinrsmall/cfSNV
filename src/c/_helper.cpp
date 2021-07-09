//
// Created by Colin Small on 7/1/21.
//

#include <vector>
#include <string>
#include <cmath>


int factorial(int n){
    if (n == 0)
        return 1;
    return n * factorial(n - 1);
}


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
std::vector<int> getIndicesWhereEqualFromString(std::string str, char comparator){
    std::vector<int> returnVec = {};

    for(int i = 0; i < str.size(); i++){
        if(str[i] == comparator)
            returnVec.push_back(i);
    }

    return returnVec;
}

double sumVector(const std::vector<double>& vec){
    double result = 0.0;
    for(auto& i : vec){
        result += i;
    }
    return result;
}

std::string toUpper(std::string s){
    for (char& c : s){
        c = toupper(c);
    }
    return s;
}

std::vector<std::string> split(const std::string& s, const std::string& delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

double meanVector(const std::vector<double>& vec){
    return sumVector(vec) / (double)vec.size();
}