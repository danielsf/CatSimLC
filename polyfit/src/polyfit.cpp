/*
Perform adapted version of the polyfit algorithm from Appendix A of
Prsa et al. (2008)
ApJ 687, 542
*/

#include "containers.h"
#define pi 3.141592654

struct Ran{

    //this structure will be based on the Xorshift random number generator
    // discovered by George Marsaglia and published in
    //Journal of Statistical Software, volume 8, no. 14 pp 1-6

    //parameters are drawn from the table on page 347 of
    //Numerical Recipes (3rd edition)
    //William H. press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
    //Cambridge University Press, 2007

    unsigned long long x;
    Ran(unsigned long long seed){

        x=seed^88172645463325252LL;
        x^=(x<<21);
        x^=(x>>35);
        x^=(x<<4);
        //printf("starting rand with %ld from seed %d\n",x,seed);
    }

    void thework(){
        x^=(x<<21);
        x^=(x>>35);
        x^=(x<<4);
    }

    double doub(){
        thework();
        return x*5.42101086242752217e-20;
    }

    int int32(){
        thework();
        int ans=int(x);
        if(ans<0)ans=-1*ans;
        return ans;
    }

};

double normal_deviate(Ran *chaos, double mu, double sig){

 int i;
 double x1,x2,y1,y2;

 x1=chaos->doub();
 x2=chaos->doub();

 y1=sqrt(-2.0*log(x1))*cos(2.0*pi*x2);
 y2=sqrt(-2.0*log(x1))*sin(2.0*pi*x2);


 return mu+y1*sig;

}


inline double power(double arg,int raised){

  //return arg raised to the integer power 'raised'

  int n;
  double ans;

  if(raised==0)return 1.0;
  else{ans=1.0;
  for(n=0;n<raised;n++){
    ans=ans*arg;
  }

  return ans;
  }

}


void naive_gaussian_solver(
array_2d<double> &aa_in, array_1d<double> &bb_in,
array_1d<double> &xx, int params){


    xx.reset_preserving_room();

    array_1d<double> buffer,bb;
    buffer.set_dim(params);
    bb.set_dim(params);

    array_2d<double> aa;
    aa.set_dim(params,params);

    buffer.set_name("naive_buffer");
    aa.set_name("naive_aa");
    bb.set_name("naive_bb");

    array_1d<int> dexes;
    dexes.set_dim(params);

    dexes.set_name("naive_dexes");

    int i,j,k;
    for(i=0;i<params;i++){
        for(j=0;j<params;j++){
            aa.set(i,j,aa_in.get_data(i,j));
        }

    }
    for(i=0;i<params;i++){
        bb.set(i,bb_in.get_data(i));
        dexes.set(i,i);
    }

    double amax,nn;
    int imax,ii;
    int row,col,rowmax,colmax;

    for(ii=0;ii<params;ii++){
        for(row=ii;row<params;row++){
	    for(col=ii;col<params;col++){
	        nn=fabs(aa.get_data(row,col));
		if((row==ii && col==ii) || nn>amax){
		
		    amax=nn;
		    rowmax=row;
		    colmax=col;
		
		}
	    }
	}
	
	if(rowmax!=ii){
	    for(i=0;i<params;i++)buffer.set(i,aa.get_data(ii,i));
	    for(i=0;i<params;i++)aa.set(ii,i,aa.get_data(rowmax,i));
	    for(i=0;i<params;i++)aa.set(rowmax,i,buffer.get_data(i));
	
	    nn=bb.get_data(ii);
	    bb.set(ii,bb.get_data(rowmax));
	    bb.set(rowmax,nn);
	
	}
	
	if(colmax!=ii){
	    for(i=0;i<params;i++)buffer.set(i,aa.get_data(i,ii));
	    for(i=0;i<params;i++)aa.set(i,ii,aa.get_data(i,colmax));
	    for(i=0;i<params;i++)aa.set(i,colmax,buffer.get_data(i));
	
	    j=dexes.get_data(ii);
	    dexes.set(ii,dexes.get_data(colmax));
	    dexes.set(colmax,j);
	
	
	}
	
	for(row=ii+1;row<params;row++){
	    nn=aa.get_data(row,ii)/aa.get_data(ii,ii);

	    for(col=0;col<params;col++){
	        aa.subtract_val(row,col,aa.get_data(ii,col)*nn);
		
	    }
	
	    bb.subtract_val(row,bb.get_data(ii)*nn);	
	}
	
	/*printf("\n");
	for(i=0;i<params;i++){
	    for(j=0;j<params;j++){
	        printf("%.4e ",aa[i*params+j]);
	    }
	    printf("\n");
	}*/
	
    }

    double err,maxerr,mindiag=-1.0,minfail,mm;
    int ifail;


    maxerr=-1.0;
    for(row=0;row<params;row++){
        for(col=0;col<params;col++){
	    if(row>col){
	        err=fabs(aa.get_data(row,col));
		if(err>maxerr)maxerr=err;
	    }
	    else if(row==col){
	        err=fabs(aa.get_data(row,col));
		if(mindiag<0.0 || err<mindiag)mindiag=err;
	    }
	}
    }


    /*if(maxerr>1.0e-6 || isnan(maxerr)){
        //printf("tridiagonalization: maxerr %e mindiag %e\n",maxerr,mindiag);
	//exit(1);
    }*/

    for(ii=params-1;ii>=0;ii--){
        buffer.set(ii,bb.get_data(ii));
	for(row=params-1;row>ii;row--){
	    buffer.subtract_val(ii,buffer.get_data(row)*aa.get_data(ii,row));
	}
	mm=buffer.get_data(ii)/aa.get_data(ii,ii);
	buffer.set(ii,mm);

    }

    for(i=0;i<params;i++){
        xx.set(dexes.get_data(i),buffer.get_data(i));
    }


    for(ii=0;ii<params;ii++){
        nn=0.0;
	for(col=0;col<params;col++){
	    nn+=xx.get_data(col)*aa_in.get_data(ii,col);
	}
	
	err=fabs(nn-bb_in.get_data(ii));
	if(bb_in.get_data(ii)!=0.0)err=err/fabs(bb_in.get_data(ii));
	if(err>maxerr || ii==0){
	    maxerr=err;
	    //if(maxerr>1.0e-6)printf("maxerr %e -- %e %e\n",maxerr,nn,bb_in.get_data(ii));
	}
    }


    if(maxerr>1.0e-3 || isnan(maxerr) || isinf(maxerr)){

	
	nn=0.0;
	minfail=-10.0;
	for(i=0;i<params;i++){
	    for(j=i+1;j<params;j++){
	       nn=0.0;
	       for(k=0;k<params;k++){
	           nn+=power(aa_in.get_data(i,k)+aa_in.get_data(j,k),2);
	       }
	       if(minfail<0.0 || nn<minfail){
	           minfail=nn;
		   ifail=j;
	       }
	    }
	}
	printf("final gaussian solver validation failed %e\n",maxerr);
	throw ifail;
    }


   //printf("naive gaussian solver maxerr %e\n",maxerr);
}


void _fit_segment_master(array_1d<double> &time, array_1d<double> &flux, array_1d<double> &sigma,
                         array_1d<double> &coeffs, int order,
                         int left_constraint, double t0, double f0,
                         int right_constraint, double t1, double f1){


    array_1d<double> bb;
    array_2d<double> mm;
    mm.set_name("mm_fit");
    bb.set_name("bb_fit");
    coeffs.reset_preserving_room();

    int i,j,k;
    mm.set_dim(order,order);
    mm.zero();
    bb.set_dim(order);
    bb.zero();

    array_1d<double> cache;
    cache.set_name("cache_fit");
    array_1d<int> cache_calc;
    cache_calc.set_name("cache_fit");
    cache.set_dim(order*2);
    cache_calc.set_dim(order*2);

    double xx,ss;
    for(k=0;k<time.get_dim();k++){
        cache_calc.zero();
        ss = power(sigma.get_data(k),2);
        for(i=0;i<order;i++){
            for(j=i;j<order;j++){
                if(cache_calc.get_data(i+j)==1){
                    xx=cache.get_data(i+j);
                }
                else{
                    if(i+j>=cache.get_dim()){
                        printf("WARNING asked for value outside of cache dim %d >= %d\n",
                        i+j,cache.get_dim());
                        exit(1);
                    }
                    xx=power(time.get_data(k),i+j);
                    cache.set(i+j,xx);
                    cache_calc.set(i+j,1);
                }
                mm.add_val(i,j,xx/ss);
            }
        }
    }

    for(i=0;i<order;i++){
        for(j=0;j<i;j++){
            if(j!=i){
                mm.set(i,j,mm.get_data(j,i));
            }
        }
    }

    for(j=0;j<time.get_dim();j++){
        xx=1.0;
        ss=power(sigma.get_data(j),2);
        for(i=0;i<order;i++){
            bb.add_val(i,flux.get_data(j)*xx/ss);
            xx*=time.get_data(j);
        }
    }

    if(left_constraint==1){
        xx=1.0;
        for(i=0;i<order;i++){
            mm.set(0,i,xx);
            xx*=t0;
        }
        bb.set(0,f0);
    }

    if(right_constraint==1){
        xx=1.0;
        for(i=0;i<order;i++){
            mm.set(1,i,xx);
            xx*=t1;
        }
        bb.set(1,f1);
    }

    naive_gaussian_solver(mm,bb,coeffs,order);

}

void _fit_segment_no_constraint(array_1d<double> &time, array_1d<double> &flux, array_1d<double> &sigma,
                                array_1d<double> &coeffs, int order){

    _fit_segment_master(time, flux, sigma,
                        coeffs, order,
                        0, -1.0, -1.0,
                        0, -1.0, -1.0);
}

void _fit_segment_one_constraint(array_1d<double> &time, array_1d<double> &flux, array_1d<double> &sigma,
                                 array_1d<double> &coeffs, int order, double t0, double f0){

    _fit_segment_master(time, flux, sigma,
                        coeffs, order,
                        1, t0, f0,
                        0, -1.0, -1.0);
}

void _fit_segment_two_constraints(array_1d<double> &time, array_1d<double> &flux, array_1d<double> &sigma,
                                  array_1d<double> &coeffs, int order,
                                  double t0, double f0, double t1, double f1){

    _fit_segment_master(time, flux, sigma,
                        coeffs, order,
                        1, t0, f0,
                        1, t1, f1);
}


double _fit_segment(array_1d<double> &time, array_1d<double> &flux, array_1d<double> &sigma,
                    array_1d<double> &t_constraint, array_1d<double> &f_constraint,
                    array_1d<double> &coeffs){
    /*
    will populate coeffs
    */


    int max_order = 11;
    int min_order = 4;
    double bic,chisq,bic_best;
    int order;

    array_1d<double> coeff_buffer;
    coeff_buffer.set_name("coeff_buffer");
    array_1d<double> f_buffer;
    f_buffer.set_name("f_buffer");
    f_buffer.set_dim(time.get_dim());

    int i,j;
    double xx;
    double ln_n = log(time.get_dim());
    double chisq_best;

    for(order=min_order;order<max_order;order++){
        if(t_constraint.get_dim()==0){
            _fit_segment_no_constraint(time,flux,sigma,coeff_buffer,order);
        }
        else if(t_constraint.get_dim()==1){
            _fit_segment_one_constraint(time,flux,sigma,coeff_buffer,order,
                                        t_constraint.get_data(0),
                                        f_constraint.get_data(0));
        }
        else if(t_constraint.get_dim()==2){
            _fit_segment_two_constraints(time,flux,sigma,coeff_buffer,order,
                                         t_constraint.get_data(0),f_constraint.get_data(0),
                                         t_constraint.get_data(1),f_constraint.get_data(1));
        }
        else{
            printf("do not know how to handle %d constraints\n",t_constraint.get_dim());
            exit(1);
        }

        f_buffer.zero();
        for(i=0;i<time.get_dim();i++){
            xx=1.0;
            for(j=0;j<order;j++){
                f_buffer.add_val(i,xx*coeff_buffer.get_data(j));
                xx*=time.get_data(i);
            }
        }
        chisq=0.0;
        for(i=0;i<time.get_dim();i++){
            chisq+=power((flux.get_data(i)-f_buffer.get_data(i))/sigma.get_data(i),2);
        }
        bic=order*ln_n+chisq;
        if(order==min_order || bic<bic_best){
            bic_best=bic;
            chisq_best=chisq;
            coeffs.reset_preserving_room();
            for(i=0;i<order;i++){
                coeffs.set(i,coeff_buffer.get_data(i));
            }
        }
    }

    /*FILE *test;
    double ff;
    if(t_constraint.get_dim()==1){
        test=fopen("junk1.txt","w");
        fprintf(test,"#order %d\n",coeffs.get_dim());
        for(i=0;i<time.get_dim();i++){
            ff=0.0;
            xx=1.0;
            for(j=0;j<coeffs.get_dim();j++){
                ff+=xx*coeffs.get_data(j);
                xx*=time.get_data(i);
            }
            fprintf(test,"%e %e %e\n",time.get_data(i),flux.get_data(i),ff);
        }
        fclose(test);
        exit(1);
    }
    if(t_constraint.get_dim()==0){
        test=fopen("junk0.txt","w");
        fprintf(test,"#order %d\n",coeffs.get_dim());
        for(i=0;i<time.get_dim();i++){
            ff=0.0;
            xx=1.0;
            for(j=0;j<coeffs.get_dim();j++){
                ff+=xx*coeffs.get_data(j);
                xx*=time.get_data(i);
            }
            fprintf(test,"%e %e %e\n",time.get_data(i),flux.get_data(i),ff);
        }
        fclose(test);
    }*/

    return chisq_best;
}


void run_polyfit(array_1d<double> &time, array_1d<double> &flux, array_1d<double> &sigma,
                 int n_knots, array_1d<double> &knots, asymm_array_2d<double> &coeffs){


    array_1d<double> seg_coeff;
    seg_coeff.set_name("seg_coeff");
    asymm_array_2d<double> coeff_buffer;
    coeff_buffer.set_name("run_polyfit_coeff_buffer");

    array_1d<double> knot_buffer;
    knot_buffer.set_name("knot_buffer");
    array_1d<double> t_buffer,f_buffer,s_buffer;
    t_buffer.set_name("t_buffer");
    f_buffer.set_name("f_buffer");
    s_buffer.set_name("s_buffer");

    array_1d<double> t_constraint, f_constraint;
    t_constraint.set_name("run_polyfit_t_constraint");
    f_constraint.set_name("run_polyfit_f_constraint");

    int n_iterations = 1000;
    double delta_knot;
    int i_iteration;
    int i;
    for(i=0;i<n_knots-1;i++){
        knot_buffer.set(i,time.get_data((i+1)*(time.get_dim()/n_knots)));
    }

    delta_knot=0.01*((time.get_data(time.get_dim()-1)-time.get_data(0))/n_knots);

    Ran dice(9123);

    int i_knot;
    int i_time_last;
    double chisq;
    double chisq_best=2.0e30;
    double xx,ff;
    printf("knot dim %d\n",knot_buffer.get_dim());
    for(i_iteration=0;i_iteration<n_iterations;i_iteration++){
        chisq = 0.0;
        coeff_buffer.reset_preserving_room();
        i_time_last=0;
        t_constraint.reset_preserving_room();
        f_constraint.reset_preserving_room();

        for(i_knot=0;i_knot<knot_buffer.get_dim();i_knot++){
            //printf("testing knot %d\n",i_knot);
            t_buffer.reset_preserving_room();
            f_buffer.reset_preserving_room();
            s_buffer.reset_preserving_room();
            for(i=i_time_last;i<time.get_dim();i++){
                if(time.get_data(i)>knot_buffer.get_data(i_knot)){
                    i_time_last=i;
                    break;
                }
                t_buffer.add(time.get_data(i));
                f_buffer.add(flux.get_data(i));
                s_buffer.add(sigma.get_data(i));
            }
            if(t_buffer.get_dim()<10){
                chisq=2.0e30;
                break;
            }

            //printf("   time steps for knot %d = %d\n",i_knot,t_buffer.get_dim());

            if(i_knot>0){
                t_constraint.set(0,knot_buffer.get_data(i_knot-1));
                ff=0.0;
                xx=1.0;
                for(i=0;i<coeff_buffer(i_knot-1)->get_dim();i++){
                    ff+=coeff_buffer.get_data(i_knot-1,i)*xx;
                    xx*=t_constraint.get_data(0);
                }
                f_constraint.set(0,ff);
            }

            if(i_knot==knot_buffer.get_dim()-1){
                t_constraint.set(1,knot_buffer.get_data(0)+time.get_data(time.get_dim()-1)-time.get_data(0));
                ff=0.0;
                xx=1.0;
                for(i=0;i<coeff_buffer(0)->get_dim();i++){
                    ff+=coeff_buffer.get_data(0,i)*xx;
                    xx*=time.get_data(0);
                }
                f_constraint.set(1,ff);
            }

            xx = _fit_segment(t_buffer, f_buffer, s_buffer,
                              t_constraint, f_constraint,
                              seg_coeff);
            chisq += xx;
            coeff_buffer.add_row(seg_coeff);
            //printf("    chisq %e %d\n",chisq,seg_coeff.get_dim());
        }

        if(i_iteration==0 || chisq<chisq_best){
            chisq_best = chisq;
            printf("set chisq_best %e per dof %e\n",chisq_best,chisq_best/time.get_dim());
            coeffs.reset_preserving_room();
            for(i=0;i<coeff_buffer.get_rows();i++){
                coeffs.add_row(coeff_buffer(i)[0]);
            }
            knots.reset_preserving_room();
            for(i=0;i<knot_buffer.get_dim();i++){
                knots.set(i,knot_buffer.get_data(i));
            }
        }

        for(i_knot=0;i_knot<knot_buffer.get_dim();i_knot++){
            xx=normal_deviate(&dice,0.0,delta_knot);
            if(i_knot==0){
                if(knot_buffer.get_data(i_knot)+xx>time.get_data(0)){
                    knot_buffer.add_val(i_knot,xx);
                }
                else{
                    i_knot--;
                }
            }
            else{
                if(knot_buffer.get_data(i_knot)+xx>knot_buffer.get_data(i_knot-1)){
                    knot_buffer.add_val(i_knot,xx);
                }
                else{
                    i_knot--;
                }
            }
        }

    }

}


int main(int iargc, char *argv[]){

    FILE *stitch_file;
    stitch_file = fopen("../workspace/validation_170824/stitched/kplr006029130_lc_stitched.txt", "r");

    array_1d<double> time,flux,sigma;
    time.set_name("time");
    flux.set_name("flux");
    sigma.set_name("sigma");

    double tt,ff,ss;
    while(fscanf(stitch_file, "%le %le %le\n", &tt,&ff,&ss)>0){
        time.add(tt);
        flux.add(ff);
        sigma.add(ss);
    }
    printf("%d %d %d\n",time.get_dim(),flux.get_dim(),sigma.get_dim());
    fclose(stitch_file);

    array_1d<double> knots;
    knots.set_name("main_knots");
    asymm_array_2d<double> coeffs;
    coeffs.set_name("main_coeffs");

    run_polyfit(time,flux,sigma,300,knots,coeffs);

}
