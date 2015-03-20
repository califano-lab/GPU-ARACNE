#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
int main()
{
    std::vector<int> nums{1,2,3};
    BOOST_FOREACH (int n, nums){
        std::cout << n << std::endl;
    }
    return 0;
}
