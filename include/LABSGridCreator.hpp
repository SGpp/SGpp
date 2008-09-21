#ifndef LABSGRIDCREATOR_HPP_
#define LABSGRIDCREATOR_HPP_


#include <iostream>
#include <sstream>

template<typename GIT, typename STORAGE>
class LABSGridCreator
{

public:
    LABSGridCreator(STORAGE &storage) : storage(storage)
    {
    }

    void createRegular(int maxlevel)
    {
        GIT index;
        for(int d = 0; d < GIT::dim; d++)
        {
            index.push(d,0,0);
        }

        //createRegularRec(GIT::dim - 1, index, 0, maxlevel - 1);
        createBorderRec(GIT::dim - 1, index, 0, maxlevel - 1);
    }

protected:

	void createBorderRec(int actDim, GIT &index, int actLevelSum, int maxLevelSum)
	{
		if(actDim == 0)
		{
			createRegular1D(index, actLevelSum, maxLevelSum);
		}
		else
		{
				typename GIT::index_type srcIndex;
            	typename GIT::level_type srcLevel;

				bool resetToZero = false;

            	index.get(actDim, srcLevel, srcIndex);

				if(srcLevel == 0)
				{
					index.push(actDim, 0, 0);
					createBorderRec(actDim-1, index, actLevelSum, maxLevelSum);
					
					index.push(actDim, 0, 1);
					createBorderRec(actDim-1, index, actLevelSum, maxLevelSum);
					
					srcLevel = 1;
					srcIndex = 1;
					index.push(actDim, srcLevel, srcIndex);
					resetToZero = true;
				}
				

				if(actLevelSum <= maxLevelSum)
				{
					createBorderRec(actDim-1, index, actLevelSum, maxLevelSum);
				}           	
				
				if(actLevelSum < maxLevelSum)
				{
					index.push(actDim, srcLevel+1, srcIndex*2 - 1);
		            createBorderRec(actDim, index, actLevelSum+1, maxLevelSum);

					index.push(actDim, srcLevel+1, srcIndex*2 + 1);
		            createBorderRec(actDim, index, actLevelSum+1, maxLevelSum);
				}
				
				if(resetToZero)
				{
					index.push(actDim, 0, 0);
				}
				else
				{
					index.push(actDim, srcLevel, srcIndex);
				}
		}
		
	}

	


	void createRegular1D(GIT &index, int actLevelSum, int maxLevelSum)
    {
        index.push(0, 0, 0);
        storage.insert(index);

        index.push(0, 0, 1);
        storage.insert(index);

        for(int l = 1; l <= maxLevelSum-actLevelSum+1; l++)
        {
            for(int i = 1; i <= 1<<(l-1); i++)
            {
                index.push(0, l, 2*i-1);
                storage.insert(index);
            }
        }

    }

private:
    STORAGE &storage;

};


#endif
