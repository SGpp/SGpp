#ifndef LINEARADAPTIVEGRIDCREATOR_HPP_
#define LINEARADAPTIVEGRIDCREATOR_HPP_

template<class GIT, class STORAGE, class ALLOC>
class LinearAdaptiveGridCreator
{

public:
    LinearAdaptiveGridCreator() : seqNumber(0)
    {

    }

    void setStorage(STORAGE* storage)
    {
        this->storage = storage;
    }


    void createRegular(int maxlevel)
    {
        GIT index;
        for(int d = 0; d < GIT::dim; d++)
        {
            index.push(d,1,1);
        }

        createRegularRec(GIT::dim - 1, index, GIT::dim, GIT::dim + maxlevel - 1);
    }

protected:
    
    void createRegular1D(GIT &index, int actLevelSum, int maxLevelSum)
    {

        for(int l = 1; l <= maxLevelSum-actLevelSum+1; l++)
        {
            for(int i = 1; i <= 1<<(l-1); i++)
            {
                GIT* insert = alloc.allocate(1);
                *insert = index;
                insert->push(0, l, 2*i-1);
                insert->rehash();
                (*storage)[insert] = nextFreeSeqNumber();
            }
        }

    }

    void createRegularRec(int actDim, GIT &index, int actLevelSum, int maxLevelSum)
    {
        if(actDim == 0 || actLevelSum == maxLevelSum)
        {
            // 1-D
            createRegular1D(index, actLevelSum, maxLevelSum);
        }
        else
        {
            createRegularRec(actDim-1, index, actLevelSum, maxLevelSum);

            typename GIT::index_type srcIndex;
            typename GIT::level_type srcLevel;

            index.get(actDim, srcLevel, srcIndex);

            index.push(actDim, srcLevel+1, srcIndex*2 - 1);
            createRegularRec(actDim, index, actLevelSum+1, maxLevelSum);

            index.push(actDim, srcLevel+1, srcIndex*2 + 1);
            createRegularRec(actDim, index, actLevelSum+1, maxLevelSum);

            index.push(actDim, srcLevel, srcIndex);

        }

    }

    unsigned int nextFreeSeqNumber()
    {
        return seqNumber++;
    }


private:
    ALLOC alloc;
    STORAGE* storage;

    unsigned int seqNumber;

};



#endif

