'use strict';

module.exports = {

    countUniqueValues: function(value) {

        if (value === null || value === undefined)
            return 0;

        const values = String(value)
            .split(',')
            .map(function(item) {
                return Number(item.trim());
            })
            .filter(function(item) {
                return Number.isFinite(item);
            });

        return new Set(values).size;
    },

    updatePlotCV: function(ui) {

    const nM = this.countUniqueValues(ui.M.value());
    const available = nM > 1;

    ui.plotCV.setPropertyValue('enable', available);

    if (!available && ui.plotCV.value())
        ui.plotCV.setValue(false);
      
    },

    view_loaded: function(ui, event) {
        this.updatePlotCV(ui);
    },

    view_updated: function(ui, event) {
        this.updatePlotCV(ui);
    },

    M_changed: function(ui, event) {

    this.updatePlotCV(ui);

    const nM = this.countUniqueValues(ui.M.value());

    if (nM > 1 && !ui.plotCV.value())
        ui.plotCV.setValue(true);
}

};