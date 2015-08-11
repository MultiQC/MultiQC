$(function () {
    $('#fq_screen_plot').highcharts({

        chart: {
            type: 'column'
        },

        title: {
            text: 'FastQ Screen Results'
        },

        xAxis: {
            categories: ['Spruce','Human','Human_chrX','Human_chrY','Mouse','Ecoli','PhiX','No Hit']
        },

        yAxis: {
            allowDecimals: false,
            min: 0,
            title: {
                text: 'Percentage Aligned'
            }
        },

        tooltip: {
            formatter: function () {
                return '<b>' + this.series.stackKey.replace('column','') + ' - ' + this.x + '</b><br/>' +
                    this.series.name + ': ' + this.y + '%<br/>' +
                    'Total: ' + this.point.stackTotal + '%';
            },
        },

        plotOptions: {
            column: {
                stacking: 'normal'
            }
        },

        series: [{
            name: 'one_hit_one_library',
            data: [0, 41.44, 0, 0, 0.01, 0, 0, 47.85],
            color: 'blue',
            stack: 'P1690_101'
        }, {
            name: 'one_hit_one_library',
            linkedTo: ':previous',
            color: 'blue',
            data: [0,41.24,0,0,0.01,0,0, 48.75],
            stack: 'P1690_102'
        }, {
            name: 'multiple_hits_one_library',
            color: '#00007f',
            data: [0,7.7,0,0,0,0,0,0],
            stack: 'P1690_101'
        }, {
            name: 'multiple_hits_one_library',
            color: '#00007f',
            linkedTo: ':previous',
            data: [0,7.29,0,0,0,0,0,0],
            stack: 'P1690_102'
        }, {
            name: 'one_hit_multiple_libraries',
            color: 'red',
            data: [0,0.89,1.1,0.02,1.08,0,0,0],
            stack: 'P1690_101'
        }, {
            name: 'one_hit_multiple_libraries',
            color: 'red',
            linkedTo: ':previous',
            data: [0,0.94,1.1,0.02,0.89,0,0,0],
            stack: 'P1690_102'
        }, {
            name: 'multiple_hits_multiple_libraries',
            color: '#7f0000',
            data: [0,2.11,0.04,0.02,0.94,0,010],
            stack: 'P1690_101'
        },{
            name: 'multiple_hits_multiple_libraries',
            color: '#7f0000',
            linkedTo: ':previous',
            data: [0,1.77,0.05,0.02,0.82,0,0,0],
            stack: 'P1690_102'
        }]
    });
});
