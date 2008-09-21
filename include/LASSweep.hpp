#ifndef LASSWEEP_HPP_
#define LASSWEEP_HPP_

template<typename GIT, typename STORAGE, typename DATA, typename FUNC>
class LASSweep
{
public:
	LASSweep(STORAGE &storage, FUNC &func) : storage(storage), func(func)
	{
	}

	LASSweep(STORAGE &storage) : storage(storage), func()
	{
	}

	void sweep(DATA &alpha)
	{
		for(int i = 0; i < GIT::dim; i++)
		{
			sweep1D(i, alpha, alpha);
		}	
	}
	

	void sweep1D(int dim, DATA &interm, DATA &alpha)
	{
		GIT index;
		
        for(int d = 0; d < GIT::dim; d++)
        {
            index.push(d,1,1);
        }
        index.rehash();
        
        if(GIT::dim == 1)
        {
    		func(storage, index, 0, interm, alpha);
        }
        else
        {
			if(dim != 0)
			{
				sweep(0, index, dim, interm, alpha);
			} 
			else
			{
				sweep(1, index, dim, interm, alpha);
			}
        }
	}

protected:
	void sweep(int recD, GIT &index, int d, DATA &interm, DATA &alpha)
	{
		typename GIT::index_type srcIndex;
        typename GIT::level_type srcLevel;

		int nextD = recD + 1;
		if(nextD == d)
		{
			nextD++;
		}

		index.get(recD, srcLevel, srcIndex);

		index.set(recD, srcLevel+1, srcIndex*2 - 1);
		if(storage.has_key(&index))
		{
			sweep(recD, index, d, interm, alpha);
		}
		
		index.set(recD, srcLevel+1, srcIndex*2 + 1);
		if(storage.has_key(&index))
		{
			sweep(recD, index, d, interm, alpha);
		}
		
		index.set(recD, srcLevel, srcIndex);
		sweepRec(nextD, index, d, interm, alpha);
	}

	void sweepRec(int recD, GIT &index, int d, DATA &interm, DATA &alpha)
	{
		if(recD >= GIT::dim)
		{
			func(storage, index, d, interm, alpha);
			index.set(d, 1, 1);
		}
		else
		{
			sweep(recD, index, d, interm, alpha);
		}
	}
	
	
	
private:
	STORAGE &storage;
	FUNC func;	
};


#endif /*LASSWEEP_HPP_*/
