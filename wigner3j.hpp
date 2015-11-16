#ifndef WIGNER3J_WIGNER3J_HPP
#define WIGNER3J_WIGNER3J_HPP

#include <cassert>
#include <cmath>
#include <vector>
#include <limits>

template<typename real_scalar_type>
int wigner3j(
    std::vector<real_scalar_type> const & thrcof,
    long const l2,
    long const l3,
    long const m2,
    long const m3,
    long const l1_min,
    long const l1_max,
    long const ndim
)
{
    int ier(0);

    // HUGE is the square root of one twentieth of the largest floating
    // point number, approximately.
    const real_scalar_type huge = std::sqrt(std::numeric_limits<real_scalar_type>::max()/real_scalar_type(20));
    const real_scalar_type srhuge = std::sqrt(huge);
    const real_scalar_type tiny = std::numeric_limits<real_scalar_type>::min();
    const real_scalar_type srtiny = std::sqrt(tiny);
    const long one(1);
    const long two(2);

    const real_scalar_type m1 = -(m2+m3);

    // Check error conditions 1 and 2
    if ( l2 < std::abs(m2) or l3 < std::abs(m3) )
    {
        std::cout<<"L2-ABS(M2) or L3-ABS(M3) 'less than zero."<<std::endl;
        ier = 1;
        return ier;
    }

    // Limits for L1
    const long l1min = std::max(std::abs(l2-l3),std::abs(m1));
    const long l1max = l2+l3;

    // If l1min=l1max, we have an analytical formula.
    if (l1min == l1max)
    {
        thrcof[0] = std::pow(real_scalar_type(-1),(real_scalar_type)std::abs(l2+m2-l3+m3) )/std::sqrt( real_scalar_type(l1min+l2+l3+one) );
        ier = 0;
        return ier;
    }
    else
    {
        long nfin = l1max - l1min + one;
        if(ndim < nfin)
        {
            // the dimension of thrcof is not large enough to hold
            // all the allowed values l1
            ier = -2;
            return;
        }

        // starting forward recursion from l1min taking nstep1 steps
        long l1 = l1min;
        real_scalar_type newfac(0);
        real_scalar_type c1(0);
        thrcof[1] = srtiny;
        real_scalar_type sum1 = real_scalar_type( two*l1 + one )*tiny;

        long lstep(1);
        while(lstep <= long(3))
        {
            lstep = lstep + long(1) ;
            l1 = l1 + long(1) ;
            real_scalar_type oldfac = newfac;
            long a1 = (l1+l2+l3 + one)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3 + one);
            long a2 = (l1+m1)*(l1-m1);
            newfac = std::sqrt( real_scalar_type(a1*a2) );

            if (l1 <= one)
            {
                //  if l1=1, (l1-1) has to be factored out of dv, so
                c1 = -real_scalar_type( (two * l1 - one)*l1*(m3-m2) )/newfac;
            }
            else
            {
                long dv = -l2*(l2 + one)*m1 + l3*(l3 + one)*m1 + l1*(l1 - one)*(m3-m2);
                real_scalar_type denom = real_scalar_type(l1 - one)*newfac;

                if (lstep > long(2) )
                {
                    real_scalar_type c1old = std::abs(c1);
                }
                else
                {
                    c1 = -real_scalar_type( (two*l1 - one)*dv )/denom;
                }
            }

            // if l1=l1min+1 the third term in the recursion eqn vanishes, hence
            real_scalar_type x = srtiny*c1;
            thrcof[1] = x;
            sum1 = sum1+tiny*real_scalar_type(two*l1 + one)*c1*c1;
            if(lstep == nfin)
            {
                real_scalar_type sumuni = sum1;
                break;
            }
        }

        real_scalar_type c2 = -real_scalar_type(l1)*oldfac/denom;

        // recursion to the next 3j-coeff x
        real_scalar_type x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)


    }

    return ier;
}
#endif //WIGNER3J_WIGNER3J_HPP
