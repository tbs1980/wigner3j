#include <iostream>
#include <wigner3j.hpp>

template<typename real_scalar_type>
void test_wigner3j()
{
    real_scalar_type l2 = 999;
    real_scalar_type l3 = 1221;
    real_scalar_type m2 = 899;
    real_scalar_type m3 = -1179;
    real_scalar_type l1 = 325;
    real_scalar_type m1 = -(m2+m3);
    std::cout<<"m1 = "<<m1<<std::endl;
    real_scalar_type l1min = std::max(std::abs(l2-l3),std::abs(m1));
    std::cout<<"l1min = "<<l1min<<std::endl;
    real_scalar_type l1max = l2 + l3;
    std::cout<<"l1max = "<<l1max<<std::endl;
    const int ndim = l1max - l1min + 1;
    std::cout<<"ndim = "<<ndim<<std::endl;
    std::cout<<"\n"<<std::endl;

}

int main(void)
{
    test_wigner3j<double>();
    return 0;
}
