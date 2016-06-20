// To serialize form into a single JSON
$.fn.serializeForm = function() {
    var o = {};
    var a = this.serializeArray();
    $.each(a, function() {
        // Take only non-empty values
        if (this.value !== "") {
            o[this.name] = this.value;
        }
    });
    return o;
};

// searchlight object
$.searchlight = new Object();
// Current query id
$.searchlight.query_id = null;
// Searchlight URL
$.searchlight.SL_URL = "http://hades:5000";
// Last displayed waveform
$.searchlight.last_wave = null;
// Waveform selection
$.searchlight.last_select = null;

// Resets state depending on the query status
$.searchlight.sl_in_query = function sl_in_query(state) {
    $("#start_query").prop("disabled", state);
    $("#cancel_query").prop("disabled", !state);
    // When starting a new query reset result URL and query ID
    if (state) {
        $.searchlight.query_id = null;
    }
};

// Create JSON for a SL request
$.searchlight.make_sl_request = function make_sl_request(query_str) {
    return JSON.stringify({query: query_str, query_id: $.searchlight.query_id});
};

// Cancel query request
$.searchlight.sl_cancel_query = function sl_cancel_query() {
    $.ajax({
        url: $.searchlight.SL_URL + "/cancel",
        type: "POST",
        data: $.searchlight.make_sl_request(null),
        contentType: "application/json; charset=utf-8",
        dataType: "json",
        cache: false,
        success: function (data) {
            console.log("Cancel request passed");
        },
        error: function (xhr, status, errorThrown) {
            console.log("Error cancelling query: " + xhr.responseText);
        }
    });
};

// Next result request
$.searchlight.sl_next_result = function sl_next_result() {
    $.ajax({
        url: $.searchlight.SL_URL + "/result",
        type: "POST",
        data: $.searchlight.make_sl_request(null),
        contentType: "application/json; charset=utf-8",
        dataType: "json",
        cache: false,
        success: function (data) {
            // check for eof
            if (data["eof"]) {
                alert("Query completed!");
                $.searchlight.sl_in_query(false);
            } else {
                var result = data["result"];
                $("#records").DataTable().row.add([result["sid"], result["id"], result["pretty_time"],
                                                  result["time"], result["len"]]).draw();
                sl_next_result();
            }
        },
        error: function (xhr, status, errorThrown) {
            alert("Error getting result: " + xhr.responseText);
            $.searchlight.sl_in_query(false);
        }
    });
};

// Request waveform (clicked from the table)
$.searchlight.sl_get_waveform = function sl_get_waveform(params) {
    $.ajax({
        url: $.searchlight.SL_URL + "/waveform",
        type: "POST",
        data: $.searchlight.make_sl_request(params),
        contentType: "application/json; charset=utf-8",
        dataType: "json",
        cache: false,
        success: function (data) {
            $.searchlight.sl_display_waveform(data["waveform"],
                "id=" + params["id"] + " time=" + params["time"]);
        },
        error: function (xhr, status, errorThrown) {
            alert("Error getting waveform: " + xhr.responseText);
        }
    });
};

// Display waveform via Flot
$.searchlight.sl_display_waveform = function sl_display_waveform(
        data, plot_label) {
    var plot_data = [];
    // Add position ticks: 0,1,2,...
    for (var i = 0; i < data.length; i++) {
        plot_data.push([i, data[i]]);
    }
    $.searchlight.last_wave = plot_data;
    $.plot("#waveform-chart", [{label: plot_label, data: plot_data}], {
        yaxis: {show: false},
        grid: {borderColor: '#ccc'},
        selection: {mode: "x"}
    });
};

// $(document).ready
$(function() {
    // Prepare table (hide tick and record columns, which are used only for waveform retrieval)
    var records_table = $("#records").DataTable({
        "columnDefs": [
        {
            "targets": [4],
            "render": function(data, type, row) {
                return data / 125; // Convert from ticks to seconds
            }
        },
        {
            "targets": [1],
            "visible": false,
            "searchable": false
        },
        {
            "targets": [3],
            "visible": false,
            "searchable": false
        }
        ]
    });

    // Toggle neighborhood elements
    $("#neighborhood").click(function() {
        var checked = this.checked;
        $(".neighb").each(function() {
            $(this).prop("disabled", !checked);
        });
    });

    // Cancel query button
    $("#cancel_query").click(function() {
        $.searchlight.sl_cancel_query();
    });

    // Table result click --- draw the waveform
    $("#records tbody").on("click", "tr", function() {
        var record = records_table.row(this).data();
        var waveform_params = {"signal": $("#signal").val(), "id": record[1],
                               "time": record[3], "len": record[4]};
        $.searchlight.sl_get_waveform(waveform_params);
    });

    // Waveform select event
    $("#waveform-chart").on("plotselected", function (event, ranges) {
        var li = Math.floor(ranges.xaxis.from);
        var ri = Math.ceil(ranges.xaxis.to);
        $.searchlight.last_select = [li, ri];
        alert("Selected: from=" + li + ", to=" + ri);
    });

    // Waveform unselect
    $("#waveform-chart").on("plotunselected", function (event) {
        $.searchlight.last_select = null;
    });

    // Submit: create JSON and send to server
    $("form").submit(function(event) {
        var form_data_json = $("form").serializeForm();
        // Disable neighborhood querying if not checked
        if (!$("#neighborhood").prop("checked")) {
            form_data_json["mimic.neighborhood.l_size"] = 0;
            form_data_json["mimic.neighborhood.r_size"] = 0;
        }
        var form_data_json_str = JSON.stringify({query: form_data_json});
        // Switch buttons
        $.searchlight.sl_in_query(true);
        // Clear the table for new results
        records_table.clear().draw();
        $.ajax({
            url: $.searchlight.SL_URL + "/query",
            type: "POST",
            data: form_data_json_str,
            dataType: "json",
            contentType: "application/json; charset=utf-8",
            cache: false,
            success: function (data) {
                // take the "tuples" field
                //data = JSON.parse(data["tuples"]);
                // Immediately query for result
                $.searchlight.query_id = data["query_id"];
                if (!$.searchlight.query_id) {
                    alert("Wrong response from the server!");
                    $.searchlight.sl_in_query(false);
                } else {
                    $.searchlight.sl_next_result();
                }
            },
            error: function (xhr, status, error) {
                alert(xhr.responseText);
                $.searchlight.sl_in_query(false);
            }
        });
        event.preventDefault();
    });
});
