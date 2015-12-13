#include <vector>

class pwlInterp2
{
    //priviate
    private:
        std::vector<double> xd, yd, zd;

        int find_posizition(std::vector<double> xxd, double xi )
        {
            int b;
            int l;
            int m;
            int r;

            if ( xi < xxd[0] || xxd[xxd.size()-1] < xi ) b = -1;
            else
            {
                l = 0;
                r = xxd.size() - 1;

                while ( l + 1 < r )
                {
                    m = ( l + r ) / 2;
                    if ( xi < xxd[m] )
                    {
                        r = m;
                    }
                    else
                    {
                        l = m;
                    }
                }
                b = l;
            }
            return b;
        }

    public:
        pwlInterp2(std::vector<double>& xd_, std::vector<double>& yd_,  std::vector<double>& zd_)
            :xd(xd_), yd(yd_), zd(zd_) {};

        std::vector<double> interp(std::vector<double>& xi_, std::vector<double>& yi_)
        {
            std::vector<double> zi(xi_.size());

            double alpha;
            double beta;
            double det;
            double dxa;
            double dxb;
            double dxi;
            double dya;
            double dyb;
            double dyi;
            double gamma;
            int i;
            int j;
            int k;

            for ( k = 0; k < zi.size(); k++ )
            {
                i = find_posizition(xd, xi_[k] );
                if ( i == -1 )
                {
                    zi[k] = 0.0;
                    continue;
                }

                j = find_posizition(yd, yi_[k] );
                if ( j == -1 )
                {
                    zi[k] = 0.0;
                    continue;
                }

                if (yi_[k] < yd[j+1] + ( yd[j] - yd[j+1] ) * (xi_[i] - xd[i] ) / ( xd[i+1] - xd[i] ) )
                {
                    dxa = xd[i+1] - xd[i];
                    dya = yd[j]   - yd[j];

                    dxb = xd[i]   - xd[i];
                    dyb = yd[j+1] - yd[j];

                    dxi =xi_[k]   - xd[i];
                    dyi =yi_[k]   - yd[j];

                    det = dxa * dyb - dya * dxb;

                    alpha = ( dxi * dyb - dyi * dxb ) / det;
                    beta =  ( dxa * dyi - dya * dxi ) / det;
                    gamma = 1.0 - alpha - beta;

                    zi[k] = alpha * zd[i+1+j*xd.size()] + beta * zd[i+(j+1)*xd.size()] + gamma * zd[i+j*xd.size()];
                }
                else
                {
                    dxa = xd[i]   - xd[i+1];
                    dya = yd[j+1] - yd[j+1];

                    dxb = xd[i+1] - xd[i+1];
                    dyb = yd[j]   - yd[j+1];

                    dxi =xi_[k]   - xd[i+1];
                    dyi =yi_[k]   - yd[j+1];

                    det = dxa * dyb - dya * dxb;

                    alpha = ( dxi * dyb - dyi * dxb ) / det;
                    beta =  ( dxa * dyi - dya * dxi ) / det;
                    gamma = 1.0 - alpha - beta;

                    zi[k] = alpha * zd[i+(j+1)*xd.size()] + beta * zd[i+1+j*xd.size()] + gamma * zd[i+1+(j+1)*xd.size()];
                }
            }
            return zi;
        }
};
