'use strict';

module.exports = {

    updateReestimationPlots: function(ui) {

    const available =
        Number(ui.twostepsN.value()) > 0;

    const controls = [
        ui.plot6,
        ui.plot7
    ];

    controls.forEach(function(control) {

        control.setPropertyValue('enable', available);

        if (!available && control.value())
            control.setValue(false);
    });
    },
    
    twostepsN_changed: function(ui, event) {
    this.updateReestimationPlots(ui);
    },

    updateSignalSupportControls: function(ui) {

    const singleSupport =
        Number(ui.supportSignalVectorN.value()) === 1;

    ui.supportSignal.setPropertyValue(
        'enable',
        singleSupport
    );

    ui.supportSignalVectorMin.setPropertyValue(
        'enable',
        !singleSupport
    );

    ui.supportSignalVectorMax.setPropertyValue(
        'enable',
        !singleSupport
    );
      
    },
    
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

    updateSupportGridPlots: function(ui) {

    const available =
        Number(ui.supportSignalVectorN.value()) > 1;

    const controls = [
        ui.plot2,
        ui.plot3,
        ui.plot4,
        ui.plot5
    ];

    controls.forEach(function(control) {

        control.setPropertyValue('enable', available);

        if (!available && control.value())
            control.setValue(false);
    });
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
        this.updateSignalSupportControls(ui);
        this.updateSupportGridPlots(ui);
        this.updateReestimationPlots(ui);
    },
    
    view_updated: function(ui, event) {
        this.updatePlotCV(ui);
        this.updateSignalSupportControls(ui);
        this.updateSupportGridPlots(ui);
        this.updateReestimationPlots(ui);
    },
    
    supportSignalVectorN_changed: function(ui, event) {
        this.updateSignalSupportControls(ui);
        this.updateSupportGridPlots(ui);
    },

    M_changed: function(ui, event) {

    this.updatePlotCV(ui);

    const nM = this.countUniqueValues(ui.M.value());

    if (nM > 1 && !ui.plotCV.value())
        ui.plotCV.setValue(true);
}

};