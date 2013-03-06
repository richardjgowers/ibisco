#! /bin/bash

echo $1 > directory

if [ ! $(ls -lh directory | awk '{ print $5 }') -eq 1 ];then

    cat directory Makefile > make_prova
        
    awk '{
            i=i+1
            if(i == 1){
                dir=$0  
            }
            regex="-L"
            where = match($0,regex)
            if(where != 0){
                print "LIBS=","-L"dir,$3, $4, $5
                    }
            else{
                if(i != 1){
                print $0
                    }
            }
        }' make_prova > Make2

    mv Make2  Makefile
    rm make_prova

fi

    rm directory 2> /dev/null

make

cp IBIsCO ~/bin/. 2> err

DIR=~/bin/.

if [ -s err ];then
    echo 'ERROR!!!'
    echo 'File IBIsCO not copied'
    cat err
else
    echo 'IBIsCO succesfull copied in directory' $DIR
fi
