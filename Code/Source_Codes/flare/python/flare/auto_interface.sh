#!/bin/bash



# generate .txt file from the comments in the header of a Fortran procedure
#
# global variables:
#
#    path - path to source file
#
# arguments:
#
#    1) procedure type (subroutine or function)
#
#    2) source file (without suffix)
#
#    3) name of the procedure & output file
#
#    4) optional prefix for the procedure name
#
autodoc() {
    procedure=$1
    src=$2
    name=$3
    prefix=$4

    # grep documentation from $src within pattern:
    # ${procedure} ${name}(
    # ...
    #   use
    doc=$(sed -n "/${procedure} ${prefix}${name}(/,\${p;/  use /q}" $src | sed '1d' | sed '$d' | sed 's/^  !//')



    cat > ${name}.txt << EOF
$doc
EOF
}



autodoc_function() { autodoc "function" $@; }
autodoc_subroutine() { autodoc "subroutine" $@; }
