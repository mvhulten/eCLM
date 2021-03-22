"""
Datm namelist and stream files
"""
from typing import List
from .namelist import Namelist, namelist_group , namelist_item

class datm_in(Namelist):

    @namelist_group
    class datm_nml():

        @namelist_item
        def anomaly_forcing(self) -> str:
            """
            If set, include anomaly forcing streams in namelist.
            """
            pass

        @namelist_item
        def bias_correct(self) -> str:
            """
            If set, include bias correction streams in namelist.
            """
            pass

        @namelist_item
        def decomp(self) -> str:
            """
            Set the decomposition option for the data model.  valid options are
            placing the global array on the root task or a simple stride-one
            load balanced one-dimensional decomposition.  other decompositions
            may be added in the future. valid values are ['root','1d'].
            1d = Vector decomposition, root = run only on the master task
            """
            pass

        @namelist_item
        def factorfn(self) -> str:
            """
            filename containing correction factors for use only with CORE2 modes (CORE2_IAF and CORE2_NYF).
            """
            pass


        @namelist_item
        def force_prognostic_true(self) -> bool:
            """
            If true, prognostic is forced to true. (default=false)
            """
            pass

        @namelist_item
        def iradsw(self) -> int:
            """
            Frequency to update radiation in number of steps (or hours if negative)
            irdasw is the radiation setting used to compute the next shortwave
            Julian date.  values greater than 1 set the next radiation to the
            present time plus 2 timesteps every iradsw.  values less than 0 turn
            set the next radiation to the present time plus two timesteps every
            -iradsw hours.  if iradsw is zero, the next radiation time is the
            present time plus 1 timestep.  (default=0.)
            """
            pass

        @namelist_item
        def presaero(self) -> bool:
            """
            If true, prescribed aerosols are sent from datm (must be true for running with CLM).
            """
            pass

        @namelist_item
        def restfilm(self) -> str:
            """
            Model restart filename for the data atmosphere model data.  This is
            optional.  If both restfils and restfilm are undefined, the restart
            filename will be read from the DATM restart pointer file (or files for multiple instances).
            """
            pass

        @namelist_item
        def restfils(self) -> str:
            """
            stream restart filename for the data atmosphere stream data.  This is
            optional.  If both restfils and restfilm are undefined, the restart
            filename will be read from the DATM restart pointer file (or files for multiple instances).
            """
            pass

        @namelist_item
        def wiso_datm(self) -> bool:
            """
            flag which if true, says to turn on water isotopes
            """
            pass

    @namelist_group
    class shr_strdata_nml():

        @namelist_item
        def datamode(self) -> str:
            """
            general method that operates on the data.  this is generally
            implemented in the data models but is set in the strdata method for
            convenience.  valid options are dependent on the data model and will
            be described elsewhere.  NULL is always a valid option and means no
            data will be generated.  default='NULL'

            datamode = "NULL"
            turns off the data model as a provider of data to the coupler.
            The atm_present flag will be set to false
            and the coupler will assume no exchange of data to or from the data model.

            datamode = "COPYALL"
            The default science mode of the data model is the COPYALL mode.
            COPYALL mode will examine the fields found in all input data streams,
            if any input field names match the field names used internally, they
            are copied into the export array and passed directly to the coupler
            without any special user code.  Any required fields not found on an
            input stream will be set to zero except for aerosol deposition fields
            which will be set to a special value.  There are several other
            scientific modes supported by the model, they are listed below.  The
            mode is selected by a character string set in the strdata namelist
            variable dataMode.

            datamode = "CORE2_NYF"
            Coordinated Ocean-ice Reference Experiments (CORE) Version 2 Normal Year Forcing.
      
            datamode = "CORE2_IAF"
            In conjunction with with CORE Version 2 atmospheric forcing data,
            provides the atmosphere forcing favored by the Ocean Model Working
            Group when coupling an active ocean model with observed atmospheric
            forcing.  This mode and associated data sets implement the CORE-IAF
            Version 2 forcing data, as developed by Large and Yeager (2008) at
            NCAR.  See the documentation for CORE version 2 datasets at
            http://data1.gfdl.noaa.gov/nomads/forms/mom4/COREv2.html Also see
            W.G.Large, S.G.Yeager (2008), The global climatology of an
            interannually varying air-sea flux data set.
            Clm Dyn doi 10.1007/s00382-008-0441-3.

            datamode = "CORE_IAF_JRA*"
            JRA55 intra-annual year forcing

            datamode = "CLMNCEP"
            In conjunction with NCEP climatological atmosphere data, provides the
            atmosphere forcing favored by the Land Model Working Group when
            coupling an active land model with observed atmospheric forcing.  This
            mode replicates code previously found in CLM (circa 2005), before the
            LMWG started using the CCSM flux coupler and data models to do
            active-land-only simulations.
            """
            pass

        @namelist_item
        def domainfile(self) -> str:
            """
            model spatial gridfile associated with the strdata.  grid information will
            be read from this file and that grid will serve as the target grid
            for all input data for this strdata input.
            """
            pass

        @namelist_item
        def dtlimit(self) -> List[float]:
            """
            array (up to 30 elements) of delta time ratio limits placed on the
            time interpolation associated with the array of streams.  this real
            value causes the model to stop if the ratio of the running maximum
            delta time divided by the minimum delta time is greater than the
            dtlimit for that stream.  for instance, with daily data, the delta
            time should be exactly one day throughout the dataset and the computed
            maximum divided by minimum delta time should always be 1.0.  for
            monthly data, the delta time should be between 28 and 31 days and the
            maximum ratio should be about 1.1.  the running value of the delta
            time is computed as data is read and any wraparound or cycling is also
            included.  this input helps trap missing data or errors in cycling.
            to turn off trapping, set the value to 1.0e30 or something similar.
            """
            pass

        @namelist_item
        def fillalgo(self) -> List[str]:
            """
            array (up to 30 elements) of fill algorithms associated with the array
            of streams.  valid options are just copy (ie. no fill), special value,
            nearest neighbor, nearest neighbor in "i" direction, or nearest
            neighbor in "j" direction.
            valid values:  'copy','spval','nn','nnoni','nnonj'
            """
            pass

        @namelist_item
        def fillmask(self) -> List[str]:
            """
            plays no role is fill algorithm at the present time.
            valid values: "nomask,srcmask,dstmask,bothmask"
            """
            pass

        @namelist_item
        def fillread(self) -> List[str]:
            """
            array (up to 30 elements) of filenames associated with the array of
            streams.  specifies the weights file to read in instead of computing the
            weights on the fly for the fill operation.  if this is set, fillalgo
            and fillmask are ignored.
            """
            pass

        @namelist_item
        def fillwrite(self) -> List[str]:
            """
            array (up to 30 elements)  of filenames associated with the array of
            streams.  specifies the weights file to generate after weights are
            computed on the fly for the fill operation.  this allows a user to
            save and reuse a set of weights later.
            """
            pass

        @namelist_item
        def mapalgo(self) -> List[str]:
            """
            array (up to 30 elements) of fill algorithms associated with the array
            of streams.  valid options are copy by index, set to special value,
            nearest neighbor, nearest neighbor in "i" direction, nearest neighbor
            in "j" direction, or bilinear.
            valid values: copy,spval,nn,nnoni,nnonj,bilinear
            """
            pass

        @namelist_item
        def mapmask(self) -> List[str]:
            """
            array (up to 30 elements) of masking algorithms for mapping input data
            associated with the array of streams.  valid options are map only from
            valid src points, map only to valid destination points, ignore all
            masks, map only from valid src points to valid destination points.
            valid values: srcmask, dstmask, nomask,bothmask
            """
            pass

        @namelist_item
        def mapread(self) -> List[str]:
            """
            array (up to 30 elements) of filenames associated with the array of
            streams.  specifies the weights file to read instead of computing
            weights on the fly for the mapping (interpolation) operation.  if this
            is set, mapalgo and mapmask are ignored.
            """
            pass

        @namelist_item
        def mapwrite(self) -> List[str]:
            """
            array (up to 30 elements) of filenames associated with the array of
            streams.  specifies the weights file to generate after weights are
            computed on the fly for the mapping (interpolation) operation.  this
            allows a user to save and reuse a set of weights later.
            """
            pass

        @namelist_item
        def readmode(self) -> List[str]:
            """
            array (up to 30 elements) of reading mode associated with the array of
            streams.  specifies the mode of reading temporal stream dataset.
            valid options are "single" (read temporal dataset one at a time) or
            "full_file" (read all entires of temporal dataset in a given netcdf file)
            valid values: single,full_file
            """
            pass

        @namelist_item
        def streams(self) -> List[str]:
            """
            character array (up to 30 elements) of stream input files.  this
            string is actually parsed by a stream method and so the format is
            specified by the stream module.  this string consists of a
            "stream_input_filename year_align year_first year_last".  the
            stream_input_filename is a stream text input file and the format and
            options are described elsewhere.  year_align, year_first, and
            year_last provide information about the time axis of the file and how
            to relate the input time axis to the model time axis.
            """
            pass

        @namelist_item
        def taxmode(self) -> List[str]:
            """
            array of time axis modes associated with the array of streams for
            handling data outside the specified stream time axis.
            valid options are to cycle the data based on the first, last, and
            align settings associated with the stream dataset, to extend the first
            and last valid value indefinitely, or to limit the interpolated data
            to fall only between the least and greatest valid value of the time array.
            valid values: cycle,extend,limit
            """
            pass

        @namelist_item
        def tintalgo(self) -> List[str]:
            """
            array (up to 30 elements) of time interpolation options associated with the array of
            streams.
            valid values: lower,upper,nearest,linear,coszen
            lower   = Use lower time-value
            upper   = Use upper time-value
            nearest = Use the nearest time-value
            linear  = Linearly interpolate between the two time-values
            coszen  = Scale according to the cosine of the solar zenith angle (for solar)
            """
            pass

        @namelist_item
        def vectors(self) -> List[str]:
            """
            list of paired colon delimited field names that should be treated as
            vectors when carrying out spatial interpolation.  unlike other
            character arrays in this namelist, this array is completely decoupled
            from the list of streams.  this is a list of vector pairs that span
            all input streams where different fields of the vector pair could
            appear in different streams.
            for example, vectors = 'u:v','taux:tauy'.
            """
            pass