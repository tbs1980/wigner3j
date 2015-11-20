#include <iostream>
#include <iomanip>
#include <wigner3j.hpp>

template<typename real_scalar_type>
void test_wigner3j()
{
    long l2 = 2;
    long l3 = 1;
    long m2 = 0;
    long m3 = 1;
    long m1 = -(m2+m3);
    long l1min =  std::max(std::abs(l2-l3),std::abs(m1));
    long l1max = l2+l3;
    long ndim = l1max-l1min+1;

    std::cout<<"l1min = "<<l1min<<std::endl;
    std::cout<<"l1max = "<<l1max<<std::endl;
    std::cout<<"ndim = "<<ndim<<std::endl;

    std::vector<real_scalar_type> thrcof = wigner3j<real_scalar_type>(l2,l3,m2,m3);

    std::cout<<"\nresuts"<<std::endl;
    for(size_t i=0;i<thrcof.size();++i)
    {
        std::cout<<i<<"\t"<<std::setprecision(20)<<thrcof[i]<<std::endl;
    }

}

int main(void)
{
    test_wigner3j<double>();
    return 0;
}
