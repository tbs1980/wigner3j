#ifndef WIGNER3J_WIGNER3J_HPP
#define WIGNER3J_WIGNER3J_HPP

#include <cassert>
#include <cmath>
#include <vector>
#include <limits>

template<typename real_scalar_type>
std::vector<real_scalar_type> wigner3j(
    long const l2, // TODO should this be unsigned long?
    long const l3, // TODO should this be unsigned long?
    long const m2,
    long const m3
)
{
    assert(std::abs(m2) <= l2);
    assert(std::abs(m3) <= l3);

    long m1 = -(m2 + m3);
    long const l1min =  std::max(std::abs(l2-l3),std::abs(m1)); // TODO should this be unsigned long?
    long const l1max = l2 + l3; // TODO should this be unsigned long?
    long const ndim = l1max - l1min + 1; // TODO should this be unsigned long?

    std::vector<real_scalar_type> thrcof(ndim);

    // HUGE is the square root of one twentieth of the largest floating
    // point number, approximately.
    real_scalar_type const huge = std::sqrt(std::numeric_limits<real_scalar_type>::max()/real_scalar_type(20));
    real_scalar_type const srhuge = std::sqrt(huge);
    real_scalar_type const tiny = std::numeric_limits<real_scalar_type>::min();
    real_scalar_type const srtiny = std::sqrt(tiny);

    long const one(1);
    long const two(2);

    // If l1min=l1max, we have an analytical formula.
    // TODO reference required here, e.g. an equation
    if (l1min == l1max)
    {
        std::cout<<"we only have one coefficeit :)"<<std::endl;
        thrcof[0] = std::pow(real_scalar_type(-1),std::abs(l2+m2-l3+m3))/std::sqrt( real_scalar_type(l1min+l2+l3+one) );
    }
    else
    {
        //std::cout<<"\n"<<std::endl;
        // starting forward recursion from l1min taking nstep1 steps
        long l1 = l1min;
        //std::cout<<"l1 = "<<l1<<std::endl;
        real_scalar_type newfac(0);
        //std::cout<<"newfac = "<<newfac<<std::endl;
        real_scalar_type c1(0);
        //std::cout<<"c1 = "<<c1<<std::endl;
        // We start with an arbitrary value.
        thrcof[0] = srtiny;
        //std::cout<<"thrcof[0] = "<<thrcof[0]<<std::endl;
        real_scalar_type sum1 = real_scalar_type(two*l1 + one)*tiny;
        //std::cout<<"sum1 = "<<sum1<<std::endl;

        // recursion starts here
        l1 = l1 + one;
        //std::cout<<"l1 = "<<l1<<std::endl;
        real_scalar_type oldfac = newfac;
        //std::cout<<"oldfac = "<<oldfac<<std::endl;
        real_scalar_type a1( (l1+l2+l3+one)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+one) );
        //std::cout<<"a1 = "<<a1<<std::endl;
        real_scalar_type a2( (l1+m1)*(l1-m1) );
        //std::cout<<"a2 = "<<a2<<std::endl;
        newfac = std::sqrt(a1*a2);
        //std::cout<<"newfac = "<<std::setprecision(20)<<newfac<<std::endl;

        assert(l1 != long(1)); // so that we check for if (l1.lt.one+eps) go to 40 in f90 version

        real_scalar_type dv( -l2*(l2 + one)*m1 + l3*(l3 + one)*m1 + l1*(l1 - one)*(m3-m2) );
        //std::cout<<"dv = "<<dv<<std::endl;
        real_scalar_type denom = real_scalar_type(l1 - one)*newfac;
        //std::cout<<"denom = "<<std::setprecision(20)<<denom<<std::endl;

        c1 = -real_scalar_type(two*l1 - one)*dv/denom;
        //std::cout<<"c1 = "<<std::setprecision(20)<<c1<<std::endl;

        // if l1=l1min+1 the third term in the recursion eqn vanishes, hence
        real_scalar_type x = srtiny*c1;
        //std::cout<<"x = "<<std::setprecision(20)<<x<<std::endl;
        thrcof[1] = x;
        //std::cout<<"thrcof[1] = "<<std::setprecision(20)<<thrcof[1]<<std::endl;
        sum1 = sum1+tiny*(two*l1 + one)*c1*c1;
        //std::cout<<"sum1 = "<<std::setprecision(20)<<sum1<<std::endl;

        if(ndim > 2)
        {
            l1 = l1 + one;
            //std::cout<<"\nl1 = "<<l1<<std::endl;
            oldfac = newfac;
            //std::cout<<"oldfac = "<<oldfac<<std::endl;
            a1 = real_scalar_type( (l1+l2+l3+one)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+one) );
            //std::cout<<"a1 = "<<a1<<std::endl;
            a2 = real_scalar_type( (l1+m1)*(l1-m1) );
            //std::cout<<"a2 = "<<a2<<std::endl;
            newfac = std::sqrt(a1*a2);
            //std::cout<<"newfac = "<<std::setprecision(20)<<newfac<<std::endl;

            assert(l1 != long(1)); // so that we check for if (l1.lt.one+eps) go to 40 in f90 version

            dv = real_scalar_type( -l2*(l2 + one)*m1 + l3*(l3 + one)*m1 + l1*(l1 - one)*(m3-m2) );
            //std::cout<<"dv = "<<dv<<std::endl;
            denom = real_scalar_type(l1 - one)*newfac;
            //std::cout<<"denom = "<<std::setprecision(20)<<denom<<std::endl;

            real_scalar_type c1old = std::abs(c1);
            //std::cout<<"c1old = "<<std::setprecision(20)<<c1old<<std::endl;

            c1 = -real_scalar_type(two*l1 - one)*dv/denom;
            //std::cout<<"c1 = "<<std::setprecision(20)<<c1<<std::endl;

            real_scalar_type c2 = -l1*oldfac/denom;

            // recursion to the next 3j-coeff x
            x = c1*thrcof[1] + c2*thrcof[0];
            //std::cout<<"x = "<<std::setprecision(20)<<x<<std::endl;
            thrcof[2] = x;
            //std::cout<<"thrcof[2] = "<<std::setprecision(20)<<thrcof[2]<<std::endl;
            real_scalar_type sumfor = sum1;
            //std::cout<<"sumfor = "<<std::setprecision(20)<<sumfor<<std::endl;
            sum1 = sum1 + (two * l1 + one)*x*x;
            //std::cout<<"sum1 = "<<std::setprecision(20)<<sum1<<std::endl;

            // see if last unnormalised 3j-coeff exceeds srhuge
            if (std::abs(x) < srhuge)
            {
                if (c1old-std::abs(c1))
                {
                    x1 = x;
                    x2 = thrcof(lstep-1)
                }
            }
        }
    }




    return thrcof;
    /*
    int ier(0);

    // HUGE is the square root of one twentieth of the largest floating
    // point number, approximately.
    real_scalar_type const huge = std::sqrt(std::numeric_limits<real_scalar_type>::max()/real_scalar_type(20));
    real_scalar_type const srhuge = std::sqrt(huge);
    real_scalar_type const tiny = std::numeric_limits<real_scalar_type>::min();
    real_scalar_type const srtiny = std::sqrt(tiny);
    long const one(1);
    long const two(2);

    real_scalar_type const m1 = -(m2+m3);

    // Check error conditions 1 and 2
    if ( l2 < std::abs(m2) or l3 < std::abs(m3) )
    {
        std::cout<<"L2-ABS(M2) or L3-ABS(M3) 'less than zero."<<std::endl;
        ier = 1;
        return ier;
    }

    // Limits for L1
    long const l1min = std::max(std::abs(l2-l3),std::abs(m1));
    long const l1max = l2+l3;

    // If l1min=l1max, we have an analytical formula.
    if (l1min == l1max)
    {
        thrcof[0] = std::pow(real_scalar_type(-1),(real_scalar_type)std::abs(l2+m2-l3+m3) )/std::sqrt( real_scalar_type(l1min+l2+l3+one) );
        ier = 0;
        return ier;
    }*/

    /*
    else
    {
        long nfin = l1max - l1min + one;
        if(ndim < nfin)
        {
            // the dimension of thrcof is not large enough to hold
            // all the allowed values l1
            ier = -2;
            return ier;
        }

        // starting forward recursion from l1min taking nstep1 steps
        long l1 = l1min;
        real_scalar_type newfac(0);
        real_scalar_type c1(0);
        thrcof[1] = srtiny;
        real_scalar_type sum1 = real_scalar_type( two*l1 + one )*tiny;
        real_scalar_type oldfac = newfac;
        real_scalar_type denom = newfac;

        long lstep(1);
        while(lstep <= long(3))
        {
            lstep = lstep + long(1) ;
            l1 = l1 + long(1) ;

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
                denom = real_scalar_type(l1 - one)*newfac;

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
        }//end of while(lstep <= long(3))

        real_scalar_type c2 = -real_scalar_type(l1)*oldfac/denom;

        // recursion to the next 3j-coeff x
        real_scalar_type x = c1*thrcof[1] + c2*thrcof[0];
        thrcof[2] = x;
        real_scalar_type sumfor = sum1;
        sum1 = sum1 + real_scalar_type(two*l1 + one)*x*x;

        if (lstep == nfin)
        {
            x1 = x;
            x2 = thrcof(lstep-1);
            x3 = thrcof(lstep-2);
            nstep2 = nfin-lstep+3;
        }


    }

    return ier;*/
}
#endif //WIGNER3J_WIGNER3J_HPP
