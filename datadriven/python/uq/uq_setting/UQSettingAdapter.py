from pysgpp import DataVector, DataMatrix
from pysgpp.extensions.datadriven.data.DataAdapter import DataAdapter
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.uq_setting import UQSetting, UQSettingFormatter
import gzip
import re

class UQSettingAdapter(DataAdapter):

    def __init__(self, uqSetting):
        """
        Constructor
        @param filename: string file name
        """
        self.uqSetting = uqSetting

    @classmethod
    def fromFile(cls, filename):
        # read from file
        s = UQSettingFormatter().deserializeFromFile(filename)
        return cls(UQSetting.fromJson(s))

    def loadData(self, qoi='_', name="train",
                 dtype=KnowledgeTypes.SIMPLE):
        """
        Reads data set from file
        @param qoi: string quatity of interest
        @param name: String for category of data set (train or test),
        default "train"
        @param dtype: knowledge type
        @return dictionary {dtype: {t: <DataContainer>}}

        WARNING: dtype parameter not supported
        """
        data = self.uqSetting.getTimeDependentResults(qoi)

        # load result in a dictionary
        ans = {}
        ans[dtype] = {}
        numDims = self.uqSetting.getDim()
        numberOfTimesteps = len(data.itervalues().next())
        # load data for all time steps
        for t, values in data.items():
            mydata = DataMatrix(numberOfTimesteps, numDims)
            sol = DataVector(size)
            for i, (sample, res) in enumerate(values.items()):
                p = DataVector(sample.getActiveUnit())
                mydata.setRow(i, p)
                sol[i] = res

            ans[dtype][t] = DataContainer(mydata, sol, name)

        return ans

    def __gzOpen(self, filename, mode="r"):
        """
        Opens a file. If the file ends with ".gz", automatically gzip
        compression is used for the file. Returns the filedescriptor
        @param filename
        @param mode default: "r" for read only
        @return file descriptor
        """
        # gzip-file
        if re.match(".*\.gz$", filename):
            # mode set for binary data?
            if not mode[-1] == "b":
                mode += "b"
            fd = gzip.open(filename, mode)
        # non gzip-file
        else:
            fd = open(filename, mode)
        return fd
