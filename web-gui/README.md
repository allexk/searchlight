Searchlight GUI
===============

Introudction
------------

Searchlight GUI is run as a Python Flask application. By default the Flaks
server is started on 0.0.0.0 at port 5000. This can be changed by editing
the searchlight.py file (statement app.run()).

Client can connect via a common HTML browser and use the GUI from there. The
GUI uses the following browser extensions:

- jQuery
- Bootstrap for the general theme
- DataTables jQuery plugin (http://datatables.net/) to display results
- Flot jQuery plugin to plot the requested waveforms

Connection Requirements
-----------------------

The GUI application assumes the queried SciDB instance is availbale at both
ports 1239 and 8080.

Port 1239 is used for the general querying, including Searchlight queries and
the waveforms queries. This is the common SciDB port for the iQuery utility,
however, this GUI uses scidb4py native Python connection instead. This
connection assumes protobuf 2.4.1 and SciDB 14.12 protocol, which is the
standara at the present.

Port 8080 is the so called HTTP SciDB shim interface. It is used solely for the
purpose of uploading auxilliary query files to the SciDB instance. This is done
via the curl utility. The GUI also supports local directories, in case the SciDB
instance is run locally and can access files from the local machine. In the
future remote SSH file copying will be added (e.g., via paramiko) to support
remotely accessible directories at the SciDB node. In the latter case the
shim interface (and the port 80880 access) will not be required.

In case SciDB is running behind a firewall, it should be possible to forward
the ports by the usual means of SSH port forwarding. For example:

    ssh -N -L localhost:8080:<scidb server>:8080 <login@ssh server)
    ssh -N -L localhost:1239:<scidb server>:1239 <login@ssh server)

Configuration
-------------

The GUI is configured by two files. The first one is flask.cfg, which configures
the Flask application. It does not usually requires any modification, unless
somehting specific is required.

The second files is confg.ini, which configures the SciDB access by the GUI.
The parameters are:

- host -- the host SciDB is running on
- user -- SciDB _HTTP shim_ user name
- passwd_file -- _HTTP shim_ password file name
    (the password will be read from the file)
- http_port -- The HTTP shim port
- scidb_port -- The main (iQuery) SciDB port
- shared_dir -- if specified, the query files will be copied locally to this
    directory, instead of the HTTP shim. The directory must be accessible
    by SciDB!

Running the GUI
---------------

After the GUI is configured, it can be run the usual way. First the Flask server
must be started:

    python searchlight.py

Then the user can use http://localhost:5000 to connect to the GUI and issue
queries. All GUI interface parameters are self-explanatory for people familiar
with Searchlight. Results might take a while to show up, depending on the query.

The results are shown in the right-side table. The user can click on any table
row, in which case the correspongin waveform is displayed below the table.
Since this requires an additional SciDB query (a simple one), this might
take some time, depending on the connection quality.

The user can cancel the current query any time. Cancel is not an immediate
call: a request is issued to the server to stop the query. When the query
is stopped, an aler will pop-up signaling the query completion. The same alert
is shown of the query finishes successfully by itself.
