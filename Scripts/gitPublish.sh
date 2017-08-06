#!/bin/bash          


#This script merges the private cellModeller repository into the public one - run it from the public repository. Run when you are sure that you want the private development stuff to get into the outside world.



while true; do
    read -p "Are you really sure about this?" yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you have both the CM4 and CM4_private repos, and are you in the public one?" yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

git remote add CM4p ../CellModeller4_private
git fetch CM4p
git branch CM4p-master remotes/CM4p/master
git merge CM4p-master
