set -x 

# You must edit the line below so that it points to the righht folder for your installation of DLIB (obtainable from dlib.net)
DLIB=/home/rejudcu/dlib-19.4
g++ -o ../bin/countVariantsInContext countVariantsInContext.cpp dcerror.cpp getSequenceFromReference.cpp hashBases.cpp -lm
g++ -o ../bin/getBackgroundCounts getBackgroundCounts.cpp dcerror.cpp hashBases.cpp -lm
g++ -o ../bin/organiseCounts organiseCounts.cpp dcerror.cpp hashBases.cpp -lm
g++ -std=c++11 -g -o ../bin/modelMutation modelMutation.cpp dcerror.cpp hashBases.cpp glfModel.cpp modelMutationFuncs.cpp runModels.cpp -lm -I $DLIB
