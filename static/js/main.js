(function() {
    'use strict';

    const Utils = {
        formatNumber: function(num, decimals) {
            if (num === null || num === undefined) return '-';
            decimals = decimals || 2;
            return Number(num).toFixed(decimals);
        },

        formatPValue: function(p) {
            if (p === null || p === undefined) return '-';
            p = Number(p);
            if (p < 0.0001) return p.toExponential(2);
            if (p < 0.001) return p.toFixed(5);
            if (p < 0.01) return p.toFixed(4);
            return p.toFixed(3);
        },

        getRegulationClass: function(regulation) {
            if (regulation === 'Upregulated') return 'text-danger';
            if (regulation === 'Downregulated') return 'text-primary';
            return 'text-muted';
        },

        getRegulationIcon: function(regulation) {
            if (regulation === 'Upregulated') return 'bi-arrow-up';
            if (regulation === 'Downregulated') return 'bi-arrow-down';
            return 'bi-dash';
        },

        showLoading: function(element, text) {
            if (!element) return;
            text = text || 'Loading...';
            element.html('<div class="text-center py-4">' +
                '<div class="spinner-border text-primary" role="status"></div>' +
                '<p class="mt-2 text-muted">' + text + '</p></div>');
        },

        showError: function(element, message) {
            if (!element) return;
            message = message || 'An error occurred';
            element.html('<div class="alert alert-danger mb-0">' +
                '<i class="bi bi-exclamation-triangle"></i> ' + message + '</div>');
        },

        showSuccess: function(element, message) {
            if (!element) return;
            message = message || 'Operation completed successfully';
            element.html('<div class="alert alert-success mb-0">' +
                '<i class="bi bi-check-circle"></i> ' + message + '</div>');
        }
    };

    const API = {
        baseUrl: '',

        get: function(endpoint, params) {
            return $.ajax({
                url: this.baseUrl + endpoint,
                method: 'GET',
                data: params || {},
                dataType: 'json'
            });
        },

        post: function(endpoint, data) {
            return $.ajax({
                url: this.baseUrl + endpoint,
                method: 'POST',
                contentType: 'application/json',
                data: JSON.stringify(data || {}),
                dataType: 'json'
            });
        },

        postFile: function(endpoint, formData) {
            return $.ajax({
                url: this.baseUrl + endpoint,
                method: 'POST',
                data: formData,
                processData: false,
                contentType: false,
                dataType: 'json'
            });
        }
    };

    const Charts = {
        instances: {},

        createDoughnut: function(canvasId, labels, data, colors) {
            const ctx = document.getElementById(canvasId);
            if (!ctx) return null;

            if (this.instances[canvasId]) {
                this.instances[canvasId].destroy();
            }

            this.instances[canvasId] = new Chart(ctx, {
                type: 'doughnut',
                data: {
                    labels: labels,
                    datasets: [{
                        data: data,
                        backgroundColor: colors,
                        borderWidth: 2,
                        borderColor: '#fff'
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: true,
                    plugins: {
                        legend: {
                            position: 'bottom',
                            labels: {
                                padding: 15,
                                font: { size: 11 }
                            }
                        }
                    },
                    cutout: '60%'
                }
            });
            return this.instances[canvasId];
        },

        createBar: function(canvasId, labels, datasets, options) {
            const ctx = document.getElementById(canvasId);
            if (!ctx) return null;

            if (this.instances[canvasId]) {
                this.instances[canvasId].destroy();
            }

            options = options || {};

            this.instances[canvasId] = new Chart(ctx, {
                type: 'bar',
                data: {
                    labels: labels,
                    datasets: Array.isArray(datasets) ? datasets : [datasets]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: true,
                    indexAxis: options.horizontal ? 'y' : 'x',
                    scales: {
                        y: {
                            beginAtZero: true,
                            ticks: { font: { size: 10 } }
                        },
                        x: {
                            ticks: { font: { size: 9 } }
                        }
                    },
                    plugins: {
                        legend: {
                            position: 'bottom',
                            labels: {
                                font: { size: 10 },
                                padding: 12
                            }
                        }
                    }
                }
            });
            return this.instances[canvasId];
        },

        destroyAll: function() {
            for (const key in this.instances) {
                if (this.instances[key]) {
                    this.instances[key].destroy();
                }
            }
            this.instances = {};
        }
    };

    const PipelineTracker = {
        updateProgress: function() {
            const items = document.querySelectorAll('.pipeline-status li');
            if (!items.length) return;

            let completed = 0;
            items.forEach(function(item) {
                if (item.classList.contains('text-success')) {
                    completed++;
                }
            });

            const total = items.length;
            const percent = Math.round((completed / total) * 100);

            const progressBar = document.getElementById('pipelineProgress');
            if (progressBar) {
                progressBar.style.width = percent + '%';
                progressBar.setAttribute('aria-valuenow', percent);
                progressBar.textContent = percent + '%';
            }

            const badge = document.getElementById('analysis-progress');
            if (badge) {
                const stages = ['Start', 'Data Ready', 'Preprocessed', 'DE Done', 'PPI Done', 'Complete'];
                badge.textContent = stages[completed] || 'Complete';
            }
        }
    };

    const GeneSearch = {
        init: function(inputId, buttonId, resultsId) {
            this.input = $('#' + inputId);
            this.button = $('#' + buttonId);
            this.results = $('#' + resultsId);
            this.timeout = null;

            if (!this.input.length) return;

            const self = this;
            this.button.on('click', function() { self.search(); });
            this.input.on('keypress', function(e) {
                if (e.which === 13) self.search();
            });
            this.input.on('input', function() {
                clearTimeout(self.timeout);
                self.timeout = setTimeout(function() {
                    if (self.input.val().length >= 2) self.search();
                }, 500);
            });
        },

        search: function() {
            const query = this.input.val().trim();
            if (!query || query.length < 1) return;

            const self = this;
            Utils.showLoading(this.results, 'Searching genes...');

            $.get('/api/gene_search?q=' + encodeURIComponent(query), function(data) {
                if (data && data.length) {
                    let html = '<div class="table-responsive"><table class="table table-sm table-hover">';
                    html += '<thead><tr><th>Gene</th><th>Log2 FC</th><th>Adj. P-value</th><th>Regulation</th></tr></thead><tbody>';
                    data.forEach(function(g) {
                        const cls = Utils.getRegulationClass(g.regulation);
                        const icon = Utils.getRegulationIcon(g.regulation);
                        html += '<tr>';
                        html += '<td><strong>' + g.gene + '</strong></td>';
                        html += '<td>' + Utils.formatNumber(g.log2_fc, 3) + '</td>';
                        html += '<td>' + g.p_adj + '</td>';
                        html += '<td class="' + cls + '"><i class="bi ' + icon + '"></i> ' + g.regulation + '</td>';
                        html += '</tr>';
                    });
                    html += '</tbody></table></div>';
                    html += '<small class="text-muted">Showing ' + data.length + ' match(s)</small>';
                    self.results.html(html);
                } else {
                    self.results.html('<div class="alert alert-info py-2"><i class="bi bi-info-circle"></i> No matching genes found</div>');
                }
            }).fail(function() {
                Utils.showError(self.results, 'Search failed. Please try again.');
            });
        }
    };

    $(document).ready(function() {
        PipelineTracker.updateProgress();

        if (document.getElementById('geneSearch') && document.getElementById('searchBtn')) {
            GeneSearch.init('geneSearch', 'searchBtn', 'searchResults');
        }

        $(document).on('click', '[data-toggle="tooltip"]', function() {
            const tooltip = new bootstrap.Tooltip(this);
            tooltip.show();
        });

        const resizeObserver = new ResizeObserver(function() {
            for (const key in Charts.instances) {
                if (Charts.instances[key]) {
                    Charts.instances[key].resize();
                }
            }
        });

        const chartsContainer = document.querySelector('.charts-container');
        if (chartsContainer) {
            resizeObserver.observe(chartsContainer);
        }
    });

    window.APP = {
        Utils: Utils,
        API: API,
        Charts: Charts,
        PipelineTracker: PipelineTracker,
        GeneSearch: GeneSearch
    };
})();
