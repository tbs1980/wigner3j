#include <iostream>
#include <wigner3j.hpp>

template<typename real_scalar_type>
void test_wigner3j()
{
    long l2 = 4;
    long l3 = 4;
    long m2 = 4;
    long m3 = 4;
    long m1 = -(m2+m3);
    long l1min =  std::max(std::abs(l2-l3),std::abs(m1));
    long l1max = l2+l3;
    long ndim = l1max-l1min+1;

    std::cout<<"l1min = "<<l1min<<std::endl;
    std::cout<<"l1max = "<<l1max<<std::endl;
    std::cout<<"ndim = "<<ndim<<std::endl;

    std::vector<real_scalar_type> thrcof = wigner3j<real_scalar_type>(l2,l3,m2,m3);

}

int main(void)
{
    test_wigner3j<double>();
    return 0;
}
