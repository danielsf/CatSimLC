/*
This file defines several classes that are used for storing vectors and
matrices.  Generally speaking, values in these functions are set using
something like

myArray.set(i,x)

to set the ith element of myArray.  They are returned using

myArray.get_data(i)

to return the ith element of myArray.  If you ask one of these classes
to return a value that doesn't exist, i.e.

myArray.set_dim(6) //setting myArray to a 0-indexed 6-element array
myArray.get_data(6) //try to return the element indexed as 6

the code will throw an exception.

Each array carries with it a member variable 'name' that is an array of char's
and can be set using

myArray.set_name("myArray")

If you choose to set 'name', then, if that array causes an exception, it will
print the name of the array to the screen, letting you know which array
was illegally manipulated.

The arrays also have a variable 'where_am_i' that can be set

myArray.set_where('inMyFunction')

This exists so that you can tell the array where it is in your code at any
given time.  'where_am_i' will also be printed to the screen in the case of an
exception.

myArray.set_where("nowhere")

sets 'where_am_i' to NULL.  Note that once 'where_am_i' is set, it persists until
it is set again, so if the array is passed between nested subroutines, it will
not necessarily know when it leaves the inner subroutine and reverts to the
outer subroutine, unless you explicitly tell it.

*/

#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define letters 500

inline int compare_char(char *s1, char *s2){
    //are two character strings the same?
    //if so, return 1
    //if not, return 0

    int i;
    //printf("\ncomparing %s %s\n",s1,s2);
    for(i=0;i<letters && (s1[i]!=0 || s2[i]!=0);i++){
     if(s1[i]!=s2[i])return 0;
    }
    return 1;
}

template <typename T>
class array_1d{

    /*
    This is a class of 1-dimensional arrays with elements of type T.

    Presently, the code will compile with T set to double or int
    To change that, change the 'template class' statements in
    containers.cpp

    */

public:

    array_1d(){
        room=0;
        dim=0;
        data=NULL;
        name=NULL;
        where_am_i=NULL;

        name_set=0;
        where_set=0;
    }

    array_1d(const array_1d<T> &in){
        room=in.room;
        dim=in.dim;
        name=NULL;
        where_am_i=NULL;
        name_set=0;
        where_set=0;
        data=new T[room];
        int i;
        for(i=0;i<room;i++){
            data[i]=in.data[i];
        }
    }

    ~array_1d(){
        if(data!=NULL){
           delete [] data;
        }

        if(name!=NULL && name_set==1){
            delete [] name;
        }

        if(where_am_i!=NULL && where_set==1){
            delete [] where_am_i;
        }

    }

    /*return a pointer to the data in this array*/
    T* get_ptr(){
        return data;
    }

    inline int contains(T val){
        int i;
        for(i=0;i<dim;i++){
            if(data[i]==val){
                return 1;
            }
        }
        return 0;
    }

    /*return the element of the array specified by int*/
    inline T get_data(int dex) const{

        if(data==NULL){
            printf("dying from get_data because data is null\n");
            die(dex);
        }

        if(dex<0 || dex>=dim){
            printf("dying from get_data because request makes no sense\n");
            die(dex);
        }

        return data[dex];

    }

    inline void set_room(int room_in){
        if(room_in<room){
            printf("dying from set_room; room is %d trying %d\n",
            room,room_in);
            die(0);
        }
        T *buffer;
        int i;

        if(data!=NULL){
            buffer=new T[dim];
            for(i=0;i<dim;i++)buffer[i]=data[i];
	        delete [] data;
	        data=new T[room_in];
	        for(i=0;i<dim;i++){
	            data[i]=buffer[i];
	    }
            delete [] buffer;
        }
        else{
            data=new T[room_in];
        }

        room=room_in;

    }

    /*add an element onto the end of the array*/
    inline void add(T in){

        if(data==NULL && dim>0){
           printf("dying from add\n");
           die(0);
        }

        if(data==NULL && room>0){
            printf("dying from add\n");
            die(0);
        }

        if(room==0 && data!=NULL){
            printf("dying from add\n");
            die(0);
        }


        if(data==NULL){
            room=10;
	    data=new T[room];
	    dim=0;
        }

        if(dim==room){
            set_room(2*room);
        }

        data[dim]=in;
        dim++;
    }

    /*set the element of the array specified by the int index to the value T.
    If you try to set an index that is beyond the current size of the array,
    zeros will be added to fill in between the current size of the array
    and the element you are setting*/
    inline void set(int dex, T val){

        int i;

        if(dex<0){
            printf("dying from set with negative dex\n");
            die(dex);
        }
        else if(dex>=dim){
            for(i=dim;i<dex+1;i++)add(0);
            set(dex,val);
        }
        else{
            data[dex]=val;
        }
    }

    /*add the value T to the element of the array indexed by int, i.e.
    array[int] = array[int] + T*/
    inline void add_val(int dex, T val){

        if(dex<0 || dex>=dim){
            printf("dying from add_val\n");
            die(dex);
        }

        data[dex]+=val;

    }

    /*subtract the value T from the element of the array indexed by int*/
    inline void subtract_val(int dex, T val){

        if(dex<0 || dex>=dim){
            printf("dying from subtract_val");
            die(dex);
        }

        data[dex]-=val;
    }

    /*divide the element of the array indexed by int by the value T, i.e.
    array[iint] = array[int]/T */
    inline void divide_val(int dex, T val){

        if(dex<0 || dex>=dim){
            printf("dying from divide_val\n");
            die(dex);
        }

        data[dex]=data[dex]/val;

    }

    /*multiply the element of the array indexed by int by the value T*/
    inline void multiply_val(int dex, T val){

        if(dex<0 || dex>=dim){
            printf("dying from multiply_val\n");
            die(dex);
        }

        data[dex]*=val;

    }


    /*set all of the elements of the array to zero*/
    void zero(){
        int i;
        for(i=0;i<room;i++)data[i]=0;
    }

    /*remove the element indexed by int from the array; shift all of the elements
    with indexes greater than int down to fill in the gap*/
    void remove(int dex){
        if(dex<0 || dex>=dim){
            return;
        }

        int i;
        for(i=dex+1;i<dim;i++){
            data[i-1]=data[i];
        }
        dim--;
    }


    /*set the length of the array to int; if the array is already carrying data
    this will not delete that data*/
    void set_dim(int ii){

        if(ii<0){
            printf("tried to set dim %d\n",dim);
            die(ii);
        }

        if(ii==0){
            reset();
            return;
        }

        T *buffer;
        int j,new_room;

        if(data!=NULL && room<ii){
            buffer=new T[room];
            for(j=0;j<room;j++)buffer[j]=data[j];
            delete [] data;
            new_room=ii;
            data=new T[new_room];
            for(j=0;j<room;j++)data[j]=buffer[j];
            delete [] buffer;
            dim=ii;
            room=new_room;
            for(;j<room;j++)data[j]=0;
        }
        else if(data!=NULL && room>=ii){
            dim=ii;
        }
        else if(data==NULL){
            data=new T[ii];
            room=ii;
            dim=ii;
            for(j=0;j<room;j++)data[j]=0;
        }

    }


    /*reduce the length of the array by one without changing the contents of the
    array*/
    void decrement_dim(){
        dim--;
        if(dim<0)dim=0;
    }

    /*increase the length of the array by one by adding a zero to the end of the
    array*/
    void increment_dim(){
        int i=dim;
        add(0);

        if(dim!=i+1){
             printf("WARNING increment_dim did not work %d %d\n",i,dim);
             die(dim);
        }
    }

    /*return the length of the array*/
    inline int get_dim() const{
       return dim;
    }

    /*set the name of the array, so that, if it causes an exception to be thrown,
    you will know which array threw the exception*/
    void set_name(char *word){
        int i,ct=0;


        if(name!=NULL && name_set==1)delete [] name;

        if(compare_char(word,"nowhere")==1){
            name_set=0;
            name=NULL;
            return;
        }

        for(i=0;word[i]!=0;i++)ct++;
        ct++;

        name=new char[ct];
        for(i=0;i<ct && word[i]!=0;i++)name[i]=word[i];
        name[ct-1]=0;

        name_set=1;

    }


    /*set where_am_i so that, if an exception is thrown, you can tell where in the program
    it happened; note that if an array is passed between multiple nested subroutines,
    this may not be the most useful message, as where_am_i does not automatically revert
    to previous settings upon leaving a subroutine*/
    void set_where(char *word) const{

        int i,ct=0;

        if(where_am_i!=NULL && where_set==1){
            delete [] where_am_i;
        }

        if(compare_char(word,"nowhere")==1){
            where_set=0;
            where_am_i=NULL;
            return;
        }

        for(i=0;word[i]!=0;i++)ct++;
        ct++;

        where_am_i=new char[ct];
        for(i=0;i<ct && word[i]!=0;i++)where_am_i[i]=word[i];
        where_am_i[ct-1]=0;

        where_set=1;

    }

    /*print the name of this array to the screen*/
    void print_name(){
        if(name!=NULL)printf("%s\n",name);
    }

    /*throw an exception; the argument is the element of the array that was asked for when
    the exception is thrown*/
    void die(int ii) const{
        printf("\nWARNING 1d array\n");

        if(name!=NULL)printf("in 1_d array %s\n",name);
        if(where_am_i!=NULL)printf("in routine %s\n",where_am_i);

        printf("asked for %d but dim %d\n",ii,dim);
        printf("room %d\n",room);

        if(data==NULL){
            printf("data is null\n\n");
        }

        int ifail=1;

        throw ifail;

    }

    /*clear the contents of the array; name and where_am_i are untouched*/
    void reset(){

        if(data==NULL && (room>0 || dim>0)){
            printf("dying from reset\n");
             die(-1);
        }

        if(data!=NULL && room==0){
            printf("dying from reset\n");
            die(-1);
        }

        delete [] data;
        data=NULL;
        room=0;
        dim=0;

    }

    void reset_preserving_room(){
        dim=0;
    }

    /*these routines exist so that array_2d can apply its name to the array_1d's that
    comprise its rows*/
    void assert_name(char *word){
        if(name_set==1){
            printf("cannot assert name; it has been set internally\n");
            die(0);
        }

        name=word;
    }

    void assert_where(char *word){
        if(where_set==1){
            printf("cannot assert where; it has been set internally\n");
            die(0);
        }

        where_am_i=word;
    }

    void assert_name_null(){
        if(name_set==1){
            delete [] name;
            name_set=0;
        }

        name=NULL;
    }

    void assert_where_null(){
        if(where_set==1){
            delete [] where_am_i;
            where_set=0;
        }

        where_am_i=NULL;
    }

    /*calculate the Euclidean norm of the array and divide all of the elements thereby.
    Return the calculated norm*/
    inline double normalize(){

        if(dim<0){
            printf("WARNING 1d array has dim %d\n",dim);
            die(-1);
        }

        if(dim==0){
            return 0.0;
        }

        double ans;
        int i;
        ans=0.0;
        for(i=0;i<dim;i++){
            ans+=data[i]*data[i];
        }

        if(ans<0.0){
            printf("WARNING square of norm %e\n",ans);

            die(-1);
        }

        if(ans>0.0){
            ans=sqrt(ans);
            for(i=0;i<dim;i++){
                data[i]=data[i]/ans;
            }
        }

        return ans;

    }

    /*return the Euclidean norm of the array without normalizing the array*/
    inline double get_norm(){

        if(dim<0){
            printf("WARNING 1d array has dim %d\n",dim);
            die(-1);
        }

        if(dim==0){
            return 0.0;
        }

        int i;
        double ans=0.0;
        for(i=0;i<dim;i++){
            ans+=data[i]*data[i];
        }
        ans=sqrt(ans);
        return ans;

    }

    /*return the square of the Euclidean norm of the array without normalizing
    the array*/
    inline double get_square_norm(){

        if(dim<0){
            printf("WARNING 1d array has dim %d\n",dim);
            die(-1);
        }

        if(dim==0){
            return 0.0;
        }

        int i;
        double ans=0.0;
        for(i=0;i<dim;i++){
            ans+=data[i]*data[i];
        }

        return ans;

    }

    /*add room for int new elements in the array*/
    void add_room(int ii){

        if(data==NULL && dim>0){
           printf("dying from add_room\n");
           die(0);
        }

        if(data==NULL && room>0){
            printf("dying from add_room\n");
            die(0);
        }

        if(room==0 && data!=NULL){
            printf("dying from add_room\n");
             die(0);
        }

        if(dim>room){
            printf("dying from add_room\n");
            die(0);
        }

        T *buffer;
        int i,old_room=room;

        if(data==NULL){
            room=ii;
            data=new T[room];
            for(i=0;i<room;i++)data[i]=0;
        }
        else{
            buffer=new T[room];
            for(i=0;i<room;i++)buffer[i]=data[i];
            delete [] data;
            room*=2;
            data=new T[room];
            for(i=0;i<old_room;i++)data[i]=buffer[i];
            delete [] buffer;
        }
}


    /*return the amount of room allotted for elements in the array.
    This is not the same as get_dim() which returns the number of occupied
    spaces in the array*/
    int get_room(){
        return room;
    }


    inline array_1d& operator=(const array_1d &in){
        if(this==&in) return *this;

        if(data!=NULL){
            delete [] data;
        }
        if(name!=NULL){
            delete [] name;
        }
        if(where_am_i!=NULL){
            delete [] where_am_i;
        }

        dim=in.dim;
        room=in.room;
        name_set=in.name_set;
        int i;
        if(in.name!=NULL){
           for(i=0;in.name[i]!=0;i++);
           name=new char[i+1];
           for(i=0;in.name[i]!=0;i++){
               name[i]=in.name[i];
           }
           name[i]=0;
        }
        if(in.where_am_i!=NULL){
            for(i=0;in.where_am_i[i]!=0;i++);
            where_am_i=new char[i+1];
            for(i=0;in.where_am_i[i]!=0;i++){
                where_am_i[i]=in.where_am_i[i];
            }
            where_am_i[i]=0;
        }

        data = new T[room];
        for(i=0;i<room;i++){
            data[i]=in.data[i];
        }
        return *this;
    }

private:

    /*this is where the data is actually stored for the array*/
    T *data;

    /*
    dim is the length of the array

    room is the length of data, which is >= dim (so that we don't have reallocate dim
    every time we add an element to the array)

    name_set is a flag that keeps track of whether name was set for this array as itself,
    or as part of an array_2d
    */
    int dim,room,name_set;

    /*where_set is like name_set for where_am_i*/
    mutable int where_set;

    /*the name of this array (for diagnostic purposes)*/
    char *name;

    /*where this array is in the program (for diagnostic purposes)*/
    mutable char *where_am_i;



};



template <typename T>
class array_2d{

public:

    /*
    Constructors.  The first constructor takes as arguments the number of
    rows and the number of columns (in that order).  You can always add rows
    to an array_2d.  Trying to add columns without first resetting the matrix
    will throw an exception.

    If you use the constructor that takes no arguments, you must specify the
    number of columns, either by using set_cols(int), set_dim(int,int) or by first adding
    a row to your blank array_2d with add_row(array_1d<T>&).
    */
    array_2d(int r,int c){
        int i,j;
        _rows=0;
        _cols=c;
        _room=r*c;

        _data=new T[_room];

        _name=NULL;
        _where_am_i=NULL;

    }

    array_2d(){
        _rows=0;
        _cols=0;
        _room=0;
        _data=NULL;
        _name=NULL;
        _where_am_i=NULL;
    }

    ~array_2d(){
        int i;

        if(_data!=NULL){
            delete [] _data;
        }

        if(_name!=NULL)delete [] _name;

        if(_where_am_i!=NULL)delete [] _where_am_i;
    }


    /*set the dimensions of the array_2d. rows first, columns second*/
    inline void set_dim(int ir, int ic){

        if(ir<0 || ic<0){
            printf("tried to set dimensions %d %d\n",ir,ic);
            die(ir,ic);
        }

        if(_data==NULL && (_rows>0 || _cols>0)){
            printf("WARNING data is null but rows %d cols %d\n",
            _rows,_cols);
            if(_name!=NULL)printf("name %s\n",_name);
            if(_where_am_i!=NULL)printf("where %s\n",_where_am_i);
            exit(1);
        }

        if(_data!=NULL && _cols<=0){
            printf("WARNING data is not null but rows %d cols %d\n",
            _rows,_cols);
            if(_name!=NULL)printf("name %s\n",_name);
            if(_where_am_i!=NULL)printf("where %s\n",_where_am_i);
            exit(1);
        }

        if(ir==_rows && ic==_cols){
            return;
        }

        if(ir==0 && ic==0){
            reset();
            return;
        }

        if((ir==0 && ic!=0) || (ic==0 && ir!=0)){
            printf("WARNING trying to set dim %d %d\n",ir,ic);
            die(ir,ic);
        }

        int i;
        if(_data!=NULL){
            delete [] _data;
        }

        _room=ir*ic;
        _rows=ir;
        _cols=ic;
        _data=new T[_room];

        int j,dex;
        for(i=0;i<_rows;i++){
            for(j=0;j<_cols;j++){
                dex=i*_cols+j;
                if(i==j){
                  _data[dex]=1;
                }
                else{
                  _data[dex]=0;
                }
            }
        }
    }

    /*set the number of columns of the array_2d. Once this is set, it cannot
    be changed without first calling reset() and deleting the contents of the
    array_2d*/
    inline void set_cols(int ii){
        reset();

        _room=100*ii;
        _rows=0;
        _cols=ii;
        _data=new T[_room];

    }


    /*return the element of the array_2d indexed by the two arguments (rows first,
    columns second)*/
    inline T get_data(int ir, int ic) const{

        if(_data==NULL){
            printf("dying from get_data\n");
            die(ir,ic);
        }
        else if(_room<_rows*_cols){
           printf("dying from get_data\n");
            die(ir,ic);
        }
        else if(ir>=_rows || ir<0 || ic>=_cols || ic<0){
           printf("dying from get_data\n");
            die(ir,ic);
        }

        return _data[ir*_cols+ic];

    }


    /*set the name member variable (for diagnostic purposes)*/
    void set_name(char *word){
        int i,ct=0;

        for(i=0;word[i]!=0;i++)ct++;
        ct++;

        if(_name!=NULL){
            delete [] _name;
        }

        _name=new char[ct];
        for(i=0;i<ct && word[i]!=0;i++)_name[i]=word[i];
        _name[ct-1]=0;

    }

    /*set where_am_i (for diagnostic purposes)*/
    void set_where(char *word) const{

        int i,ct=0;

        for(i=0;word[i]!=0;i++)ct++;
        ct++;

        if(_where_am_i!=NULL)delete [] _where_am_i;
        _where_am_i=new char[ct];

        for(i=0;i<ct && word[i]!=0;i++)_where_am_i[i]=word[i];
        _where_am_i[ct-1]=0;

    }

    /*print name to screen*/
    void print_name(){
        if(_name!=NULL)printf("%s\n",_name);
    }


    void set_row_room(int row_room_in){

        if(row_room_in*_cols<_room){
            printf("dying from set_row_room; row_roomis %d trying to set %d\n",
            _room/_cols, row_room_in);
            die(0,0);
        }

        if(_cols<=0){
           printf("dying from set_row_room because cols are %d\n",_cols);
           die(0,0);
        }

        T *buffer;
        int i,j;

        if(_data!=NULL){
            buffer=new T[_room];
            for(i=0;i<_room;i++){
                buffer[i]=_data[i];
            }
            delete[] _data;

            _data=new T[row_room_in*_cols];
            for(i=0;i<_room;i++){
               _data[i]=buffer[i];
            }
            delete [] buffer;
            _room=row_room_in*_cols;
            for(;i<_room;i++){
                _data[i]=0;
            }
        }
        else{
            _data = new T[row_room_in*_cols];
        }

        _room=row_room_in*_cols;
    }


    /*add the array_1d as a row to this array_2d.  If this array_2d is blank,
    then the number of columns in this array_2d will be set to the length
    of the input array_1d*/
    void add_row(const array_1d<T> &in){

        if(_data==NULL && _room>0){
            printf("dying from add_row\n");
            die(0,0);
        }

        if(_data==NULL && _rows>0){
            printf("dying from add_row\n");
            die(1,0);
        }

        if(_data==NULL && _cols>0){
            printf("dying from add_row\n");
            die(2,0);
        }

        if(_room<_rows*_cols){
            printf("dying from add_row\n");
            die(3,0);
        }

        if((_cols<=0 || _room<=0) && _data!=NULL){
            printf("dying from add_row\n");
            die(4,0);
        }

        if(_cols>0 && _cols!=in.get_dim()){
            printf("dying from add_row cols are: %d but got: %d\n",_cols,in.get_dim());
            die(-2,0);
        }

        if(_cols==0){
            _cols=in.get_dim();
        }

        if(_cols==0){
            printf("about add a row but cols zero\n");
            die(-1,-1);
        }

        int i,j;

        if(_data==NULL){
            set_row_room(100);
        }

        if(_rows==_room/_cols){
            set_row_room(2*_room/_cols);
        }

        for(i=0;i<_cols;i++){
            try{
                _data[_rows*_cols+i]=in.get_data(i);
            }
            catch(int iex){
                printf("dying from array_2d add_row when actually copying data\n");
                die(_rows,i);
            }
        }

        _rows++;
    }

    /*set the row indexed by int to the array_1d provided*/
    void set_row(int dex, const array_1d<T> &in){

        if(dex<0){
            printf("tried to set to negative row\n");
            die(dex,0);
        }

        if(_data==NULL && _room>0){
            printf("dying from set_row\n");
            die(0,0);
        }

        if(_data==NULL && _rows>0){
             printf("dying from set_row\n");
            die(1,0);
        }

        if(_data==NULL && _cols>0){
             printf("dying from set_row\n");
            die(2,0);
        }

        if(_room<_rows){
            printf("dying from set_row\n");
            die(3,0);
        }

        if((_cols<=0 || _room<=0) && _data!=NULL){
             printf("dying from set_row\n");
            die(4,0);
        }

        if(_cols>0 && _cols!=in.get_dim()){
             printf("dying from set_row\n");
            printf("columns do not match\n");
            die(-2,0);
        }

        if(dex<0){
             printf("dying from set_row\n");
            die(dex,0);
        }

        int i;
        if(dex>=_rows){
            for(i=_rows;i<dex+1;i++)add_row(in);
        }
        else{
            for(i=0;i<_cols;i++){
                try{
                    _data[dex*_cols+i]=in.get_data(i);
                }
                catch(int ifail){
                    die(dex,_cols);
                }
            }
        }

    }

    /*set the element indexed by the two ints to the value provided.
    If you try to set a row that is beyond the current size of this
    array_2d, the code will add rows of just zero onto the end of
    the array_2d until there is room.  If you try to set a column
    that is beyond the current size of this array_2d, the code
    will throw an exception.*/
    inline void set(int ir, int ic, T val){

        if(ir<0){
            printf("tried to set to negative row\n");
            die(ir,ic);
        }

        if(ic<0 || ic>=_cols){
            printf("dying from set\n");
            die(ir,ic);
        }

        if(_cols<=0){
            printf("\nYou cannot use set(int,int) on a 2d array if cols are zero\n");
            die(ir,ic);
        }

        if(_data==NULL){
            printf("dying from set\n");
            die(ir,ic);
        }

        int i;
        array_1d<T> vector;
        if(ir>=_rows){
            for(i=0;i<_cols;i++)vector.set(i,0);
            while(_rows<=ir)add_row(vector);

        }

        _data[ir*_cols+ic]=val;

    }

    /*set all of the elements of this array_2d to zero*/
    void zero(){
        int i;
        for(i=0;i<_room;i++){
            _data[i]=0;
        }
    }

    /*
    add the provided value to the indexed element, i.e.
    array[int1][int2] = array[int1][int2]+T
    */
    inline void add_val(int ir, int ic, T val){

        if(ir>=_rows || ic>=_cols || _data==NULL || ir<0 || ic<0){
            printf("dying from add_val\n");
             die(ir,ic);
        }
        _data[ir*_cols+ic]+=val;
    }

    /*subtract the provided value from the indexed element*/
    inline void subtract_val(int ir, int ic, T val){

        if(ir>=_rows || ic>=_cols || _data==NULL || ir<0 || ic<0){
            printf("dying from subtract_val\n");
            die(ir,ic);
        }
        _data[ir*_cols+ic]-=val;
    }

    /*multiply the indexed element by the provided value*/
    inline void multiply_val(int ir, int ic, T val){

        if(ir>=_rows || ic>=_cols || _data==NULL || ir<0 || ic<0){
            printf("dying from multiply_val\n");
            die(ir,ic);
        }
        _data[ir*_cols+ic]*=val;
    }

    /*divide the indexed element by the provided value*/
    inline void divide_val(int ir, int ic, T val){

        if(ir>=_rows || ic>=_cols || _data==NULL || ir<0 || ic<0){
            printf("dying from divide_val\n");
            die(ir,ic);
        }
        _data[ir*_cols+ic]/=val;

    }

    /*reset the contents of the array_2d; name and where_am_i are untouched*/
    void reset(){

        //printf("resetting %s\n",name);

        //set_name("resetting");

        int i;

        if(_data==NULL && (_rows>0 || _cols>0 || _room>0)){
            printf("resetting but data is null and something is wrong\n");
            die(-1,-1);
        }

        if(_room==0 && _data!=NULL){
            die(-1,-1);
        }

        if(_cols==0 && _data!=NULL){
            die(-1,-1);
        }

        if(_room<_rows){
            die(-2,-2);
        }

        if(_data!=NULL){
            delete [] _data;

            _data=NULL;
            _room=0;
            _rows=0;
            _cols=0;
        }

    }


    void reset_preserving_room(){
        _rows=0;
    }

    /*reduce the number of rows by one; the contents of the array_2d
    are untouched*/
    void decrement_rows(){

        if(_rows==0){
            printf("WARNING trying to decrement rows but rows already zero\n");
            die(0,0);
        }

        _rows--;
    }

    /*return the number of rows*/
    inline int get_rows() const{
        return _rows;
    }

    /*return the number of columns*/
    inline int get_cols() const{
        return _cols;
    }

    /*throw an exception; the arguments are for indicating which element
    the code tried to access when the exception was thrown*/
    void die(int ir, int ic) const{
        printf("\nWARNING 2d array\n");

        if(_name!=NULL)printf("in 2d_array %s\n",_name);
        if(_where_am_i!=NULL)printf("in routine %s\n",_where_am_i);

        printf("asked for %d %d\n",ir,ic);
        printf("but dimensions are %d %d\n",_rows,_cols);
        printf("room %d\n",_room);

        if(_data==NULL){
            printf("data is null\n");
        }

        int ifail=1;

        throw ifail;
    }

    /*remove the row indexed by the int.  All of the rows with
    indexes greater than the argument are shifted to fill in the gap*/
    void remove_row(int dex){

        if(dex<0 || dex>=_rows)return;

        int i,j;
        for(i=dex+1;i<_rows;i++){
            for(j=0;j<_cols;j++){
                _data[(i-1)*_cols+j]= _data[i*_cols+j];
            }
        }

        _rows--;

    }

    /*return a pointer to the indexed array, i.e.

    myArray2d(i) is a pointer to the ith row of myArray2d, so that

    *myArray2d(i) behaves just like an array_1d
    */
    inline array_1d<T> operator()(int dex) const{

        if(dex<0 || dex>=_rows){
            printf("WARNING asked for row %d but only have %d\n",dex,_rows);
            die(-1,-1);
        }
        array_1d<T> vv;
        int i;
        vv.set_dim(_cols);
        for(i=0;i<_cols;i++){
            vv.set(i,_data[dex*_cols+i]);
        }
        return vv;

    }


    inline double normalize_row(int dex){
        if(dex>=_rows){
            printf("WARNING trying to normalize row %d; only have %d\n",
                   dex,_rows);

            die(dex,0);
        }
        double ans=0.0;
        int i;
        for(i=0;i<_cols;i++){
            ans+=_data[dex*_cols+i]*_data[dex*_cols+i];
        }
        ans=sqrt(ans);
        for(i=0;i<_cols;i++){
            _data[dex*_cols+i]/=ans;
        }
        return ans;
    }

private:

   int _rows,_cols,_room;

   T *_data;

   char *_name;

   mutable char *_where_am_i;

};

template <typename T>
class asymm_array_2d{

    /*
    This is a class for a 2-dimensional 'matrix' in which each row has a
    different number of columns.

    The rows are stored as array_1d's.

    As before, the code will currently compile so that you can have
    asymm_array_2d's of either ints or doubles.
    */

public:
    asymm_array_2d(){
        name=NULL;
        where_am_i=NULL;
        data=NULL;
        rows=0;
        row_room=0;
    }

    ~asymm_array_2d(){
        int i;

        for(i=0;i<row_room;i++){
            try{
                data[i].assert_name_null();
            }
            catch(int iex){
                printf("in asymm 2d destructor\n");
                die(0);
            }

            try{
                data[i].assert_where_null();
            }
            catch(int iex){
                printf("in asymm 2d destructor\n");
                die(0);
            }
        }

        //printf("calling 2d destructor on %s\n",name);

        if(data!=NULL){
            delete [] data;
        }

        if(name!=NULL)delete [] name;

        if(where_am_i!=NULL)delete [] where_am_i;


    }


    /*set the name of this asymm_array_2d for diagnostic purposes*/
    void set_name(char *word){
        int i,ct=0;

        for(i=0;word[i]!=0;i++)ct++;
        ct++;

        if(name!=NULL){
            delete [] name;
        }

        name=new char[ct];
        for(i=0;i<ct && word[i]!=0;i++)name[i]=word[i];
        name[ct-1]=0;

        for(i=0;i<row_room;i++){
            try{
                data[i].assert_name(name);
            }
            catch(int iex){
                printf("in asymm 2d set name\n");
                die(0);
            }
        }

    }


    /*set the location of this asymm_array_2d for diagnostic purposes*/
    void set_where(char *word) const{

        int i,ct=0;

        for(i=0;word[i]!=0;i++)ct++;
        ct++;

        if(where_am_i!=NULL)delete [] where_am_i;
        where_am_i=new char[ct];

        for(i=0;i<ct && word[i]!=0;i++)where_am_i[i]=word[i];
        where_am_i[ct-1]=0;

        for(i=0;i<row_room;i++){
            try{
                data[i].assert_where(where_am_i);
            }
            catch(int iex){
                printf("in asymm 2d set where\n");
                die(0);
            }
        }

    }


    /*add a row to the end of this asymm_array_2d*/
    void add_row(const array_1d<T> &in){

        if(data==NULL){
            row_room=2;
            rows=0;
            data=new array_1d<T>[row_room];
        }

        array_1d<T> *buffer;
        int i,j;

        if(rows==row_room){
            buffer=new array_1d<T>[rows];
            for(i=0;i<rows;i++){
                buffer[i].set_dim(data[i].get_dim());
                for(j=0;j<data[i].get_dim();j++){
                    buffer[i].set(j,data[i].get_data(j));
                }
            }
            delete [] data;

            i=row_room/2;
            if(i<100)i=100;

            row_room*=2;
            data=new array_1d<T>[row_room];

            for(i=0;i<rows;i++){
                data[i].set_dim(buffer[i].get_dim());
                for(j=0;j<buffer[i].get_dim();j++){
                    data[i].set(j,buffer[i].get_data(j));
                }
            }
            delete [] buffer;
        }

        data[rows].set_dim(in.get_dim());
        for(i=0;i<in.get_dim();i++){
            data[rows].set(i,in.get_data(i));
        }
        rows++;

        for(i=0;i<row_room;i++){
                try{
                    data[i].assert_name(name);
                }
                catch(int iex){
                    printf("in asymm 2d add row (asserting name)\n");
                    die(0);
                }

                try{
                    data[i].assert_where(where_am_i);
                }
                catch(int iex){
                    printf("in asymm 2d add row (asserting where)\n");
                    die(0);
                }
        }

    }


    /*set the indexed row to the provided array_1d;
    If you set a row beyond the current size of this asymm_array_2d,
    then empty rows will be used to fill in the gaps
    */
    void set_row(int ir, const array_1d<T> &vv){

        int i;
        array_1d<T> empty;

        while(rows<=ir){
            add_row(empty);
        }

        data[ir].reset();
        for(i=0;i<vv.get_dim();i++){
            data[ir].set(i,vv.get_data(i));
        }

    }

    /*
    remove the row indexed by int.  Rows with indexes greater than the provided
    index will be shifted down to fill in the gap
    */
    void remove_row(int dex){

        if(dex<0 || dex>=rows){
            printf("WARNING asking to remove %d from asymm\n",dex);
            die(dex);
        }

        int i,j;
        for(i=dex;i<rows-1;i++){
            data[i].set_dim(data[i+1].get_dim());
            for(j=0;j<data[i+1].get_dim();j++){
                data[i].set(j,data[i+1].get_data(j));
            }
        }
        data[rows-1].reset();
        rows--;

    }

    /*set all of the elements in this asymm_array_2d to zero*/
    void zero(){
        int i;
        for(i=0;i<row_room;i++)data[i].zero();
    }

    /*
    set the element indexed by the ints to the provided value.

    If you specify an element beyond the present size of the asymm_array_2d,
    zeros will be used to fill in the empty space.

    Note that there is no longer any restriction on the specified column
    index as each row in asymm_array_2d is allowed to have a different
    number of columns.
    */
    inline void set(int ir, int ic, T val){

        array_1d<T> empty;
        int i;

        while(rows<=ir){
            add_row(empty);
        }

        data[ir].set(ic,val);
    }


    /*
    Return the indexed element
    */
    inline T get_data(int ir, int ic) const{

        if(ir<0 || ir>=rows){
            printf("WARNING asking for asymm 2d data %d %d but rows %d\n",
	    ir,ic,rows);
	    die(ir);
        }


        try{
           return data[ir].get_data(ic);
        }
        catch(int iex){
            printf("tried to get asymm 2d data %d %d\n",ir,ic);
	    die(ir);
        }

        return data[ir].get_data(ic);
    }

    /*add T to the end of the row specified by int*/
    inline void add(int dex, T val){
        int i;
        if(dex<0){
            printf("in asymm 2d add\n");
            die(dex);
        }

        array_1d<T> empty;

        while(dex>=rows){
            add_row(empty);
        }
        data[dex].add(val);
    }

    /*add the value T onto the indexed element, i.e.
    asymm_array[int1][int2] = asymm_array[int1][int2] + T
    */
    inline void add_val(int ir, int ic, T val){
        if(ir<0 || ir>=rows){
            printf("in asymm 2d add_val\n");
	    die(ir);
        }

        data[ir].add_val(ic,val);
    }

    /*subtract the value T from the indexed element*/
    inline void subtract_val(int ir, int ic, T val){
        if(ir<0 || ir>=rows){
            printf("in asymm 2d subtract_val\n");
            die(ir);
        }

        data[ir].subtract_val(ic,val);
    }

    /*divide the indexed element by the value T*/
    inline void divide_val(int ir, int ic, T val){
        if(ir<0 || ir>=rows){
            printf("in asymm 2d divide_val\n");
            die(ir);
        }

        data[ir].divide_val(ic,val);
    }

    /*multiply the indexed elmement by the value T*/
    inline void multiply_val(int ir, int ic, T val){
        if(ir<0 || ir>=rows){
            printf("in asymm 2d multiply_val\n");
            die(ir);
        }

        data[ir].multiply_val(ic,val);
    }

    /*replace the indexed row with the provided array_1d*/
    void replace_row(int dex, array_1d<T> &pt){
        if(dex<0 || dex>=rows){
            printf("WARNING trying to replace row %d in asymm, but only have %d\n",
            dex,rows);

            die(dex);
        }

        data[dex].reset();
        int i;
        for(i=0;i<pt.get_dim();i++){
            data[dex].add(pt.get_data(i));
        }

    }

    inline int get_rows() const{
        return rows;
    }

    inline int get_cols(int dex) const{

        if(data==NULL){
            return 0;

            //printf("WARNING asking for cols in asymm array 2d\n");
            //die(dex);
        }

        if(dex<0){
            printf("WARNING asking for cols in asymm array 2d\n");
            die(dex);
        }

        if(dex>=rows){
            return 0;
        }

        return data[dex].get_dim();

    }

    /*throw an exception; the argument indicates the row index being
    called for when the exception was thrown*/
    void die(int ir) const{
        printf("\nWARNING asymm 2d array\n");

        if(name!=NULL)printf("in 2d_array %s\n",name);
        if(where_am_i!=NULL)printf("in routine %s\n",where_am_i);

        printf("asked for %d\n",ir);
        printf("but dimensions are %d\n",rows);
        printf("row_room %d\n",row_room);

        if(data==NULL){
            printf("data is null\n");
        }

        int ifail=1;

        throw ifail;
    }

    /*reset the contents of this asymm_array_2d*/
    void reset(){

        //printf("resetting %s\n",name);

        //set_name("resetting");

        int i;

        if(data==NULL && (rows>0 || row_room>0)){
            printf("resetting but data is null and something is wrong\n");
            die(-1);
        }

        if(row_room==0 && data!=NULL){
            die(-1);
        }

        if(row_room<rows){
            die(-2);
        }

        if(data!=NULL){
            delete [] data;
            data=NULL;
            row_room=0;
            rows=0;

        }

    }

    void reset_preserving_room(){
        int i;
        for(i=0;i<rows;i++){
            data[i].reset_preserving_room();
        }
        rows=0;
    }

    /*return a pointer to the row indexed by int, i.e.

    myAsymmArray(i) is a pointer to the array_1d in the ith row

    so

    *myAsymmArray(i) behaves like an array_1d

    */
    inline array_1d<T>* operator()(int dex){

        if(dex<0 || dex>=rows){
            printf("WARNING asked for row %d but only have %d\n",dex,rows);
        }

        return &data[dex];
    }

private:

    /*
    rows is the number of rows

    row_room is the number of array_1d's allocated for rows (so
    the code knows how much room it has when you call add_row)
    */
    int rows,row_room;

    /*this will be allocated as an array of array_1d's for storing the
    rows of this asymm_array_2d*/
    array_1d<T> *data;

    /*the name of this asymm_array_2d for diagnostic purposes*/
    char *name;

    /*the location of this asymm_array_2d in the code, for diagnostic purposes*/
    mutable char *where_am_i;

    /*name and where_am_i will be passed down to the array_1d's making up the rows
    of this assymm_array_2d*/
};


///////////////*BELOW ARE ROUTINES FOR SORTING array_1d's/////////////

/*a driver for a merge sort as described here

http://en.wikipedia.org/wiki/Merge_sort

The user should not call this routine; it is called by sort_and_check, which
the user should call
*/
template <typename T>
void merge_sort(const array_1d<T>&,array_1d<int>&,int,int);



/*
The first array_1d<T> is the input data

The second array_1d<T> will store the sorted output

The array_1d<int> are unique id's associated with each element of the input data.
It will be rearranged to preserve the relationship between these indexes and the
sorted data.

This routine will check to make sure that the output is sorted from lowest
to highest and that the unique id's are preserved relative to the input data.

The routine returns the maximum error fabs(input(unique_id=i)-output(unique_id=i)).

If this maximum error is greater than 10^-12, then the code throws an exception
*/
template <typename T>
double sort_and_check(const array_1d<T>&, array_1d<T>&, array_1d<int>&);

template <typename T>
void sort(const array_1d<T>&, array_1d<T>&, array_1d<int>&);


/*return the index of the element of the array_1d that is closest in value to T;

ASSUMES THAT THE array_1d IS SORTED FROM LOWEST TO HIGHEST
*/
template <typename T>
int get_dex(const array_1d<T>&, T);

#endif
