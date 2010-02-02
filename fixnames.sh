#!/bin/sh

git filter-branch --env-filter '

n=$GIT_AUTHOR_NAME
m=$GIT_AUTHOR_EMAIL

case ${GIT_AUTHOR_NAME} in
        amrit)    n="Amrit Jalan";        m="amrit@mit.edu" ;;
        cfgold)   n="Franklin Goldsmith"; m="cfgold@mit.edu" ;;
        gberan)   n="Gregory Beran";      m="gberan@mit.edu" ;;
        gmagoon)  n="Greg Magoon";        m="gmagoon@mit.edu" ;;
        jdmo)     n="Jeffrey Mo";         m="jdmo@mit.edu" ;;
        jwallen)  n="Josh Allen";         m="jwallen@mit.edu" ;;
        karma)    n="Karma James";        m="karma09@mit.edu" ;;
        mrharper) n="Michael Harper";     m="mrharper@mit.edu" ;;
        petway)   n="Sally Petway";       m="petway@mit.edu" ;;
        rwashcra) n="Robert Ashcraft";    m="rwashcra@mit.edu" ;;
        rwest)    n="Richard West";       m="rwest@mit.edu" ;;
        sandeeps) n="Sandeep Sharma";     m="sandeeps@mit.edu" ;;
esac

export GIT_AUTHOR_NAME="$n"
export GIT_AUTHOR_EMAIL="$m"
export GIT_COMMITTER_NAME="$n"
export GIT_COMMITTER_EMAIL="$m"
'