const getTextDimensions = function (text, fontSize) {
  let div = document.createElement("div");

  div.innerText = text;
  div.style.position = "absolute";
  div.style.float = "left";
  div.style.fontSize = fontSize;
  div.style.whiteSpace = "nowrap";
  div.style.visibility = "hidden";

  document.body.append(div);
  let width = div.clientWidth;
  let height = div.clientHeight;
  div.remove();

  return [width, height];
};

const getActiveCaller = () => {
  return +d3.select("input[name=dataset]:checked").node().value;
};

const plotChromosomeView = function (plotData) {
  let data = plotData;
  const height = 400;
  const width = 800;
  const margin = {
    top: 10,
    right: 30,
    bottom: 40,
    left: 60,
    between: 40,
  };
  const plotHeight = (height - margin.top - margin.between - margin.bottom) / 2;

  const xScale = d3
    .scaleLinear()
    .domain([0, data.length])
    .range([0, width - margin.left - margin.right]);

  const yScale = d3.scaleLinear().range([plotHeight, 0]);

  const yScaleVAF = d3.scaleLinear().domain([0, 1]).range([plotHeight, 0]);

  const xAxis = (g) => g.call(d3.axisBottom(xScale).ticks(5));

  const yAxis = (g) =>
    g.call(
      d3
        .axisLeft(yScale)
        .ticks(8)
        .tickFormat((y, i) => (i % 2 == 0 ? y : ""))
    );

  const yAxisVAF = (g) => g.call(d3.axisLeft(yScaleVAF).ticks(5));

  const svg = d3
    .select("#chromosome-view")
    .attr("preserveAspectRatio", "xMinYMin meet")
    .attr("viewBox", [0, 0, width, height])
    .attr(
      "style",
      "max-width: 100%; height: auto; max-height: 500px; height: intrinsic;"
    );

  // Clip path for annotations
  svg
    .append("clipPath")
    .attr("id", "annotation-clip")
    .append("rect")
    .attr("width", width - margin.left - margin.right)
    .attr("height", 2 * plotHeight + margin.between);

  // Clip path for log ratios
  svg
    .append("clipPath")
    .attr("id", "lr-area-clip")
    .append("rect")
    .attr("width", width - margin.left - margin.right)
    .attr("height", plotHeight);

  const plotArea = svg
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`)
    .attr("class", "plot-area");

  const lrArea = plotArea
    .append("g")
    .attr("id", "lr-plot")
    .attr("clip-path", "url(#lr-area-clip)");
  const vafArea = plotArea
    .append("g")
    .attr("id", "vaf-plot")
    .attr("clip-path", "url(#lr-area-clip)")
    .attr("transform", `translate(0, ${plotHeight + margin.between})`);

  const regions = lrArea.append("g").attr("class", "regions");

  const segments = lrArea.append("g").attr("class", "segments");

  const lrGrid = lrArea.append("g").attr("class", "grid");
  const vafGrid = vafArea.append("g").attr("class", "grid");

  const plotGrid = function () {
    lrGrid
      .selectAll(".gridline")
      .data(yScale.ticks())
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", xScale(0))
            .attr("x2", xScale(data.length))
            .attr("y1", (d) => yScale(d))
            .attr("y2", (d) => yScale(d))
            .attr("class", "gridline"),
        (update) =>
          update.call((update) =>
            update
              .transition()
              .attr("y1", (d) => yScale(d))
              .attr("y2", (d) => yScale(d))
          ),
        (exit) => exit.remove()
      );

    vafGrid
      .selectAll(".gridline")
      .data(yScaleVAF.ticks())
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", xScale(0))
            .attr("x2", xScale(data.length))
            .attr("y1", yScaleVAF)
            .attr("y2", yScaleVAF)
            .attr("class", "gridline"),
        (update) =>
          update.call((update) =>
            update.transition().attr("y1", yScaleVAF).attr("y2", yScaleVAF)
          ),
        (exit) => exit.remove()
      );
  };

  const plotRegions = function (regionData) {
    regions
      .selectAll(".point")
      .data(regionData, (d) => [d.start, d.end, d.log2])
      .join(
        (enter) =>
          enter
            .append("circle")
            .attr("class", "point")
            .attr("cx", (d) => xScale(d.start))
            .attr("cy", (d) => yScale(d.log2))
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0)
            .call((enter) =>
              enter.transition().duration(500).attr("fill-opacity", 0.3)
            ),
        (update) =>
          update.attr("fill-opacity", 0.3).call((update) =>
            update
              .transition()
              .attr("cx", (d) => xScale(d.start))
              .attr("cy", (d) => yScale(d.log2))
          ),
        (exit) => exit.transition().attr("fill-opacity", 0).remove()
      );
  };

  const plotSegments = function (segmentData) {
    segments
      .selectAll(".segment")
      .data(segmentData, (d) => [d.start, d.end])
      .join(
        (enter) =>
          enter
            .append("path")
            .attr("class", "segment")
            .attr(
              "d",
              (d) =>
                `M${xScale(d.start)} ${yScale(d.log2)} L ${xScale(
                  d.end
                )} ${yScale(d.log2)}`
            )
            .attr("stroke", "orange")
            .attr("stroke-width", 2)
            .attr("stroke-opacity", 0)
            .call((enter) => enter.transition().attr("stroke-opacity", 1)),
        (update) =>
          update
            .attr("stroke-opacity", 1)
            .call((update) =>
              update
                .transition()
                .attr(
                  "d",
                  (d) =>
                    `M${xScale(d.start)} ${yScale(d.log2)} L ${xScale(
                      d.end
                    )} ${yScale(d.log2)}`
                )
            ),
        (exit) => exit.transition().attr("stroke-opacity", 0).remove()
      );
  };

  const plotAnnotations = function () {
    plotArea
      .selectAll(".annotation")
      .attr("clip-path", "url(#annotation-clip)")
      .data(data.annotations, (d) => [data.chromosome, d.name, d.start, d.end])
      .join(
        (enter) => {
          let annotation_group = enter.append("g").attr("class", "annotation");
          annotation_group
            .append("rect")
            .attr("class", "annotation-marker")
            .attr("x", (d) => xScale(d.start))
            .attr("width", (d) => xScale(d.end) - xScale(d.start))
            .attr("height", height - margin.top - margin.bottom)
            .attr("stroke", "#000")
            .attr("stroke-width", 0.5)
            .attr("fill", "#333")
            .attr("fill-opacity", 0)
            .attr("pointer-events", "none")
            .call((enter) => enter.transition().attr("fill-opacity", 0.05));
          annotation_group
            .append("rect")
            .attr("class", "annotation-label-background")
            .attr("x", (d) => {
              let [labelWidth, labelHeight] = getTextDimensions(
                d.name,
                "0.8rem"
              );
              return (
                xScale(d.start + (d.end - d.start) / 2) - labelWidth / 2 - 5
              );
            })
            .attr("y", (d) => {
              let [labelWidth, labelHeight] = getTextDimensions(
                d.name,
                "0.8rem"
              );
              return plotHeight + margin.between / 2 - labelHeight / 2 - 2;
            })
            .attr("width", (d) => getTextDimensions(d.name, "0.8rem")[0] + 10)
            .attr("height", (d) => getTextDimensions(d.name, "0.8rem")[1] + 4)
            .attr("fill", "#EEE")
            .attr("rx", 4);
          annotation_group
            .append("text")
            .attr("class", "annotation-label")
            .text((d) => d.name)
            .attr("x", (d) => xScale(d.start + (d.end - d.start) / 2))
            .attr("y", plotHeight + margin.between / 2)
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "central");
          return annotation_group;
        },
        (update) => {
          update
            .selectAll(".annotation-label")
            .call((update) =>
              update
                .transition()
                .attr("x", (d) => xScale(d.start + (d.end - d.start) / 2))
            );
          update.selectAll(".annotation-label-background").call((update) =>
            update.transition().attr("x", (d) => {
              let [labelWidth, labelHeight] = getTextDimensions(
                d.name,
                "0.8rem"
              );
              return (
                xScale(d.start + (d.end - d.start) / 2) - labelWidth / 2 - 5
              );
            })
          );
          update
            .selectAll(".annotation-marker")
            .attr("fill-opacity", 0.05)
            .call((update) =>
              update
                .transition()
                .attr("x", (d) => xScale(d.start))
                .attr("width", (d) => xScale(d.end) - xScale(d.start))
            );
          return update;
        },
        (exit) =>
          exit
            .transition()
            .attr("fill-opacity", 0)
            .attr("stroke-opacity", 0)
            .remove()
      );
  };

  const plotVAF = function () {
    vafArea
      .selectAll(".point")
      .data(data.vaf, (d) => [data.chromosome, d.pos, d.vaf])
      .join(
        (enter) =>
          enter
            .append("circle")
            .attr("class", "point")
            .attr("cx", (d) => xScale(d.pos))
            .attr("cy", (d) => yScaleVAF(d.vaf))
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0)
            .call((enter) => enter.transition().attr("fill-opacity", 0.3)),
        (update) =>
          update.call((update) =>
            update.transition().attr("cx", (d) => xScale(d.pos))
          ),
        (exit) => exit.transition().attr("fill-opacity", 0).remove()
      );
  };

  svg
    .append("g")
    .attr("transform", `translate(${margin.left}, ${margin.top})`)
    .attr("class", "zoom-layer")
    .append("rect")
    .attr("width", width - margin.left - margin.right)
    .attr("height", height - margin.top - margin.bottom)
    .attr("fill", "transparent")
    .call(
      d3
        .drag()
        .on("start", (e) => {
          d3.select(".zoom-layer")
            .append("rect")
            .attr("class", "zoom-region")
            .attr("pointer-events", "none")
            .attr("x", e.x)
            .attr("width", 0)
            .attr("height", height - margin.bottom - margin.top)
            .attr("stroke-width", 0)
            .attr("fill", "#333")
            .attr("fill-opacity", 0.1);
        })
        .on("drag", (e) => {
          let leftBound = Math.min(e.x, e.subject.x);
          let width = Math.abs(Math.max(0, e.x) - e.subject.x);
          if (leftBound + width > xScale.range()[1]) {
            width = xScale.range()[1] - leftBound;
          }
          d3.select(".zoom-region")
            .attr("x", Math.max(0, Math.min(e.x, e.subject.x)))
            .attr("width", width);
        })
        .on("end", (e) => {
          d3.select(".zoom-region").remove();
          let xMin = Math.max(0, Math.min(e.x, e.subject.x));
          let xMax = Math.min(xScale.range()[1], Math.max(e.x, e.subject.x));
          if (xMax - xMin < 3) {
            return;
          }
          zoomToRegion(
            data.chromosome,
            xScale.invert(xMin),
            xScale.invert(xMax),
            0
          );
        })
    )
    .on("click", (e) => {
      let [xMin, xMax] = xScale.domain();
      if (xMax - xMin !== data.length) {
        // Only reset if something actually changed
        update(data, getActiveCaller());
      }
    });

  svg
    .append("g")
    .attr("transform", `translate(${margin.left},${height - margin.bottom})`)
    .attr("class", "x-axis")
    .call(xAxis);
  svg
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`)
    .attr("class", "y-axis")
    .call(yAxis);
  svg
    .append("g")
    .attr(
      "transform",
      `translate(${margin.left},${margin.top + plotHeight + margin.between})`
    )
    .attr("class", "y-axis")
    .call(yAxisVAF);

  // Labels
  svg
    .append("text")
    .attr("transform", `translate(${width / 2},${height})`)
    .attr("class", "x-label")
    .text(data.label)
    .attr("text-anchor", "middle");

  svg
    .append("text")
    .attr(
      "transform",
      `translate(0,${margin.top + plotHeight / 2}) rotate(-90)`
    )
    .attr("class", "y-label")
    .text("log2 ratio")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge");

  svg
    .append("text")
    .attr(
      "transform",
      `translate(0,${height - margin.bottom - plotHeight / 2}) rotate(-90)`
    )
    .attr("class", "y-label")
    .text("VAF")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge");

  const update = function (plotData, xMin, xMax) {
    data = plotData;
    const callerIndex = getActiveCaller();
    const fitToData = d3.select("#chromosome-fit-to-data").node().checked;
    xScale.domain([0, data.length]);
    if (xMin && xMax) {
      xScale.domain([Math.max(xMin, 0), Math.min(xMax, data.length)]);
      const yMin = data.callers[callerIndex].ratios
        .filter((d) => d.start > xMin && d.start < xMax)
        .map((d) => d.log2)
        .reduce((a, d) => (d < a ? d : a));
      const yMax = data.callers[callerIndex].ratios
        .filter((d) => d.start > xMin && d.start < xMax)
        .map((d) => d.log2)
        .reduce((a, d) => (d > a ? d : a));
      if (fitToData) {
        d3.selectAll(".data-range-warning").classed("hidden", true);
        const padding = (yMax - yMin) * 0.05;
        yScale.domain([yMin - padding, yMax + padding]);
      } else {
        d3.selectAll(".data-range-warning").classed(
          "hidden",
          !(yMin < -2 || yMax > 2)
        );
        yScale.domain([-2, 2]);
      }
    } else {
      const [yMin, yMax] = d3.extent(
        data.callers[callerIndex].ratios,
        (d) => d.log2
      );
      if (fitToData) {
        d3.selectAll(".data-range-warning").classed("hidden", true);
        const padding = (yMax - yMin) * 0.05;
        yScale.domain([yMin - padding, yMax + padding]);
      } else {
        d3.selectAll(".data-range-warning").classed(
          "hidden",
          !(yMin < -2 || yMax > 2)
        );
        yScale.domain([-2, 2]);
      }
    }
    svg.transition().select(".x-axis").duration(500).call(xAxis);
    svg.transition().select(".y-axis").duration(500).call(yAxis);
    svg.select(".x-label").text(data.label);
    plotGrid();
    plotRegions(data.callers[callerIndex].ratios);
    plotVAF();
    plotSegments(data.callers[callerIndex].segments);
    plotAnnotations();
  };

  const getZoomRange = () => {
    return xScale.domain();
  };

  update(data);

  return {
    update: update,
    getZoomRange: getZoomRange,
  };
};

const plotGenomeView = () => {
  const totalLength = d3.sum(cnvData.map((d) => d.length));
  const height = 400;
  const width = 800;
  const margin = {
    top: 10,
    right: 30,
    bottom: 60,
    left: 60,
    between: 20,
  };
  const panelWidths = cnvData.map(
    (d) => ((width - margin.left - margin.right) * d.length) / totalLength
  );
  const panelHeight =
    (height - margin.top - margin.bottom - margin.between) / 2;

  const xScales = cnvData.map((d, i) =>
    d3.scaleLinear().domain([0, d.length]).range([0, panelWidths[i]])
  );

  const yScaleRange = 2;
  const yScale = d3
    .scaleLinear()
    .domain([-yScaleRange, yScaleRange])
    .range([panelHeight, 0]);

  const yScaleVAF = d3.scaleLinear().domain([0, 1]).range([panelHeight, 0]);

  const yAxis = (g) => g.call(d3.axisLeft(yScale).ticks(5));

  const yAxisVAF = (g) => g.call(d3.axisLeft(yScaleVAF).ticks(5));

  const svg = d3
    .select("#genome-view")
    .attr("preserveAspectRatio", "xMinYMin meet")
    .attr("viewBox", [0, 0, width, height])
    .attr("style", "max-width: 100%; max-height: 500px; height: auto;");

  const plotArea = svg
    .append("g")
    .attr("transform", `translate(${margin.left}, ${margin.top})`);

  const lrArea = plotArea.append("g").attr("class", "genome-view-area");
  const vafArea = plotArea
    .append("g")
    .attr("class", "genome-view-area")
    .attr("transform", `translate(0,${panelHeight + margin.between})`);

  const addPanels = function (g) {
    const panels = g
      .selectAll(".chromosome-panel")
      .data(cnvData)
      .join("g")
      .attr("data-index", (d, i) => i)
      .attr("class", "chromosome-panel")
      .attr(
        "transform",
        (d, i) =>
          `translate(${i === 0 ? 0 : d3.sum(panelWidths.slice(0, i))}, 0)`
      );

    // Panel backgrounds
    panels
      .append("rect")
      .attr("class", "bg-rect")
      .attr("width", (d, i) => panelWidths[i])
      .attr("height", panelHeight)
      .attr("fill", "#FFF")
      .attr("stroke", "#333");

    return panels;
  };

  const lrPanels = addPanels(lrArea);
  const vafPanels = addPanels(vafArea);

  const lrGrid = lrPanels
    .append("g")
    .attr("class", "grid")
    .attr("data-index", (d, i) => i);
  const vafGrid = vafPanels
    .append("g")
    .attr("class", "grid")
    .attr("data-index", (d, i) => i);

  const ratioPanels = lrPanels
    .append("g")
    .attr("class", "regions")
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("data-index", (d, i) => i);

  const segmentPanels = lrPanels
    .append("g")
    .attr("class", "segments")
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("data-index", (d, i) => i);

  const plotRegions = () => {
    ratioPanels
      .selectAll(".point")
      // Only plot every fifth point for performance
      .data(
        (d) =>
          d.callers[getActiveCaller()].ratios.filter((r, i) => i % 5 === 0),
        (d) => d.start
      )
      .join(
        (enter) =>
          enter
            .append("circle")
            .attr("class", "point")
            .attr("cx", (d, i, g) =>
              xScales[g[i].parentNode.dataset.index](d.start)
            )
            .attr("cy", yScale(-yScaleRange - 0.2))
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0.3)
            .call((enter) =>
              enter.transition().attr("cy", (d) => yScale(d.log2))
            ),
        (update) =>
          update.call((update) =>
            update
              .transition()
              .attr("cx", (d, i, g) =>
                xScales[g[i].parentNode.dataset.index](d.start)
              )
              .attr("cy", (d) => yScale(d.log2))
          ),
        (exit) =>
          exit
            .transition()
            .attr("cy", yScale(yScaleRange + 0.2))
            .remove()
      );
  };

  const plotSegments = () => {
    segmentPanels
      .selectAll(".segment")
      // Only draw segments that will actually be visible
      .data(
        (d) =>
          d.callers[getActiveCaller()].segments.filter(
            (s) => s.end - s.start > totalLength / width
          ),
        (d) => [d.start, d.end, d.log2]
      )
      .join(
        (enter) =>
          enter
            .append("path")
            .attr("class", "segment")
            .attr("d", (d, i, g) => {
              let j = g[i].parentNode.dataset.index;
              let xScale = xScales[j];
              return `M${xScale(d.start)} ${yScale(
                -yScaleRange - 0.2
              )} L ${xScale(d.end)} ${yScale(-yScaleRange - 0.2)}`;
            })
            .attr("stroke", "orange")
            .attr("stroke-width", 2)
            .call((enter) =>
              enter.transition().attr("d", (d, i, g) => {
                let j = g[i].parentNode.dataset.index;
                let xScale = xScales[j];
                return `M${xScale(d.start)} ${yScale(d.log2)} L ${xScale(
                  d.end
                )} ${yScale(d.log2)}`;
              })
            ),
        (update) =>
          update.attr("d", (d, i, g) => {
            let j = g[i].parentNode.dataset.index;
            let xScale = xScales[j];
            return `M${xScale(d.start)} ${yScale(d.log2)} L ${xScale(
              d.end
            )} ${yScale(d.log2)}`;
          }),
        (exit) =>
          exit
            .transition()
            .attr("d", (d, i, g) => {
              let j = g[i].parentNode.dataset.index;
              let xScale = xScales[j];
              return `M${xScale(d.start)} ${yScale(
                yScaleRange + 0.2
              )} L ${xScale(d.end)} ${yScale(yScaleRange + 0.2)}`;
            })
            .remove()
      );
  };

  // Log ratio grid lines
  lrGrid
    .selectAll(".gridline")
    .data(yScale.ticks())
    .join("line")
    .attr("x1", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[0])
    .attr("x2", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[1])
    .attr("y1", (d) => yScale(d))
    .attr("y2", (d) => yScale(d))
    .attr("class", "gridline");

  // VAF grid lines
  vafGrid
    .selectAll(".gridline")
    .data(yScaleVAF.ticks())
    .join("line")
    .attr("x1", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[0])
    .attr("x2", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[1])
    .attr("y1", (d) => yScaleVAF(d))
    .attr("y2", (d) => yScaleVAF(d))
    .attr("class", "gridline");

  // VAF
  vafPanels
    .append("g")
    .attr("class", "vaf")
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("data-index", (d, i) => i)
    .selectAll(".point")
    .data((d) => d.vaf)
    .join("circle")
    .attr("cx", (d, i, g) => xScales[g[i].parentNode.dataset.index](d.pos))
    .attr("cy", (d) => yScaleVAF(d.vaf))
    .attr("r", 2)
    .attr("fill", "#333")
    .attr("fill-opacity", 0.3);

  // Clip path to create inner stroke
  const overlayClip = d3.selectAll(".genome-view-area").append("g");
  overlayClip
    .selectAll(".panel-overlay-clip")
    .data(panelWidths)
    .join("clipPath")
    .attr("class", "panel-overlay-clip")
    .attr("id", (d, i) => `panel-${i}-overlay-clip`)
    .append("rect")
    .attr("width", (d) => d)
    .attr("height", panelHeight);

  let selectedChromosomeIndex = 0;
  const selectChromosome = (chromosome) => {
    previousChromosomeIndex = selectedChromosomeIndex;
    selectedChromosomeIndex = cnvData.findIndex(
      (d) => d.chromosome === chromosome
    );
    if (previousChromosomeIndex === selectedChromosomeIndex) {
      return;
    }
    plotArea.selectAll(".panel-overlay").classed("selected", false);
    plotArea
      .selectAll(`.panel-${selectedChromosomeIndex}-overlay`)
      .classed("selected", true);
    chromosomeView.update(cnvData[selectedChromosomeIndex]);
  };

  const overlays = d3.selectAll(".genome-view-area").append("g");
  overlays
    .selectAll(".panel-overlay")
    .data(panelWidths)
    .join("rect")
    .attr(
      "class",
      (d, i) => `panel-overlay panel-${i}-overlay${i === 0 ? " selected" : ""}`
    )
    .attr(
      "transform",
      (d, i) => `translate(${i === 0 ? 0 : d3.sum(panelWidths.slice(0, i))}, 0)`
    )
    .attr("data-index", (d, i) => i)
    .attr("width", (d) => d)
    .attr("height", panelHeight)
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("fill", "#000")
    .attr("fill-opacity", 0)
    .attr("stroke", "forestgreen")
    .on("mouseenter", (e) => {
      plotArea.selectAll(".panel-overlay").attr("fill-opacity", 0);
      d3.selectAll(`.panel-${e.target.dataset.index}-overlay`).attr(
        "fill-opacity",
        0.2
      );
    })
    .on("mouseout", (e) => {
      d3.selectAll(`.panel-${e.target.dataset.index}-overlay`).attr(
        "fill-opacity",
        0
      );
    })
    .on("click", (e, d, i) =>
      selectChromosome(cnvData[e.target.dataset.index].chromosome)
    );

  // Y axes
  svg
    .append("g")
    .attr("transform", `translate(${margin.left}, ${margin.top})`)
    .attr("class", "y-axis")
    .call(yAxis);

  svg
    .append("g")
    .attr(
      "transform",
      `translate(${margin.left}, ${margin.top + panelHeight + margin.between})`
    )
    .attr("class", "y-axis")
    .call(yAxisVAF);

  // Labels
  vafPanels
    .append("text")
    .attr(
      "transform",
      (d, i) =>
        `translate(${panelWidths[i] / 2},${panelHeight + 10}) rotate(-90)`
    )
    .attr("class", "x-label")
    .text((d, i) => d.label)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "central");

  svg
    .append("text")
    .attr(
      "transform",
      `translate(0,${margin.top + panelHeight / 2}) rotate(-90)`
    )
    .attr("class", "y-label")
    .text("log2 ratio")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge");

  svg
    .append("text")
    .attr(
      "transform",
      `translate(0,${
        margin.top + margin.between + (3 * panelHeight) / 2
      }) rotate(-90)`
    )
    .attr("class", "y-label")
    .text("VAF")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge");

  const update = () => {
    plotRegions();
    plotSegments();
  };

  const getSelectedChromosome = () => selectedChromosomeIndex;

  update();

  return {
    update: update,
    getSelectedChromosome: getSelectedChromosome,
    selectChromosome: selectChromosome,
  };
};

const populateTable = () => {
  const table = d3.select("#cnv-table");
  if (table.empty()) {
    return null;
  }

  const tableHeader = table.append("thead").append("tr");
  const tableBody = table.append("tbody");
  let applyFilter = d3.select("#table-filter-toggle").node().checked;

  let tableData = null;

  const getColumnDef = (col) => {
    switch (col) {
      // Integers
      case "start":
      case "length":
      case "end":
        return {
          class: "right",
          format: (x) => x.toLocaleString(undefined, {}),
        };

      // Floating point, allow for missing numbers
      case "copyNumber":
      case "corrCopyNumber":
      case "baf":
      case "cn":
        return {
          class: "right tooltip-trigger",
          format: (x) => {
            if (x !== null && !isNaN(x)) {
              return x.toLocaleString(undefined, { minimumFractionDigits: 2 });
            }
            return "NA";
          },
        };

      case "view":
        return {
          class: "view-region-link",
          format: (x) => x,
        };

      case "type":
        return {
          class: "left",
          format: (x) => x,
        };

      // Strings
      case "caller":
      case "chromosome":
      case "gene":
      default:
        return {
          class: "left",
          format: (x) => x,
        };
    }
  };

  const getColumnLabel = (col) => {
    const columns = {
      view: "View",
      caller: "Caller",
      chromosome: "Chr",
      genes: "Genes",
      start: "Start",
      end: "End",
      length: "Length",
      type: "Type",
      cn: "CN",
      baf: "BAF",
    };

    if (columns[col]) {
      return columns[col];
    }

    return col;
  };

  const tooltip = d3
    .select(".container")
    .append("div")
    .attr("class", "copy-number-tooltip hidden")
    .call((d) =>
      d
        .append("table")
        .call((t) =>
          t
            .append("thead")
            .append("tr")
            .selectAll("th")
            .data(["caller", "type", "cn"])
            .join("th")
            .text(getColumnLabel)
        )
        .call((t) => t.append("tbody"))
    );

  const copyNumberTooltip = (data, x, y) => {
    tooltip.classed("hidden", false);
    tooltip
      .select("tbody")
      .selectAll("tr")
      .data(data)
      .join("tr")
      .selectAll("td")
      .data((d) => Object.entries(d))
      .join("td")
      .text((d) => getColumnDef(d[0]).format(d[1]));
  };

  const showCopyNumberTooltip = (e) => {
    const tableRow = e.target.parentNode.dataset.index;
    const others = tableData[tableRow].others;

    copyNumberTooltip(others);
  };

  const hideCopyNumberTooltip = (e) => {
    tooltip.classed("hidden", true);
  };

  const update = () => {
    const callerIndex = getActiveCaller();
    tableData = cnvData
      .map((d) =>
        d.callers[callerIndex].cnvs
          .filter((di) => !applyFilter || (applyFilter && di.passed_filter))
          .map((di) => {
            const allCols = { view: "🔍", chromosome: d.chromosome, ...di };
            // Don't display caller and filter status in table
            const { caller, passed_filter, ...cols } = allCols;
            return cols;
          })
      )
      .flat();

    if (tableData.length === 0) {
      tableData = [{ "No data to display": [] }];
    } else {
      // Find the corresponding copy numbers from the other caller(s)
      for (cnv of tableData) {
        // Same chromosome
        let chromData = cnvData.filter((d) => d.chromosome === cnv.chromosome);
        // ... different caller
        let callerData = chromData[0].callers.filter(
          (_d, i) => i !== callerIndex
        );
        // ... same gene
        let otherCnvs = callerData
          .map((d) =>
            d.cnvs
              .filter(
                (c) => cnv.genes.filter((g) => c.genes.includes(g)).length > 0
              )
              .map((c) => {
                return { caller: d.name, type: c.type, cn: c.cn };
              })
          )
          .flat();
        cnv.others = otherCnvs;
      }
    }

    tableHeader
      .selectAll("th")
      .data(
        Object.keys(tableData[0]).filter((k) => k !== "others"),
        (d) => d
      )
      .join("th")
      .text(getColumnLabel)
      .attr("class", (d) => getColumnDef(d).class);

    tableBody
      .selectAll("tr")
      .data(tableData)
      .join("tr")
      .attr("data-chromosome", (d) => d.chromosome)
      .attr("data-start", (d) => d.start)
      .attr("data-length", (d) => d.length)
      .attr("data-index", (_, i) => i)
      .selectAll("td")
      .data((d) => Object.entries(d).filter(([k, _]) => k !== "others"))
      .join("td")
      .text(([key, value]) => getColumnDef(key).format(value))
      .attr("class", ([key, _]) => getColumnDef(key).class);

    tableBody.selectAll(".view-region-link").on("click", (e) => {
      const dataset = e.target.parentElement.dataset;
      zoomToRegion(
        dataset.chromosome,
        Number(dataset.start),
        Number(dataset.start) + Number(dataset.length)
      );
    });

    tableBody
      .selectAll(".tooltip-trigger")
      .on("mouseenter", showCopyNumberTooltip)
      .on("mouseout", hideCopyNumberTooltip)
      .on("mousemove", (e) => {
        const dims = tooltip.node().getBoundingClientRect();
        const offset = d3.select(".container").node().getBoundingClientRect().y;
        const maxHeight =
          document.documentElement.clientHeight - offset - dims.height;
        tooltip
          .style("top", `${Math.min(e.layerY, maxHeight)}px`)
          .style("left", `${e.layerX - dims.width - 20}px`);
      });
  };

  d3.select("#table-filter-toggle").on("change", (event) => {
    applyFilter = event.target.checked;
    update();
  });

  const getData = () => tableData;

  update();

  return {
    getData: getData,
    update: update,
  };
};

// Populate dataset picker
d3.select("#dataset-picker")
  .selectAll("div")
  .data(cnvData[0].callers)
  .join("div")
  .call((e) => {
    e.append("input")
      .attr("type", "radio")
      .property("checked", (_, i) => i === 0)
      .attr("value", (_, i) => i)
      .attr("id", (d) => `dataset-${d.name}`)
      .attr("name", "dataset");
    return e;
  })
  .call((e) => {
    e.append("label")
      .attr("for", (d) => `dataset-${d.name}`)
      .text((d) => d.label);
    return e;
  });

const zoomToRegion = (chromosome, start, end, padding = 0.05) => {
  genomeView.selectChromosome(chromosome);
  let bpPadding = (end - start) * padding;
  chromosomeView.update(
    cnvData[genomeView.getSelectedChromosome()],
    start - bpPadding,
    end + bpPadding
  );
};

const chromosomeView = plotChromosomeView(cnvData[0]);
const genomeView = plotGenomeView();
const resultsTable = populateTable();

d3.select("#chromosome-fit-to-data").on("change", (e) => {
  let [xMin, xMax] = chromosomeView.getZoomRange();
  chromosomeView.update(
    cnvData[genomeView.getSelectedChromosome()],
    xMin,
    xMax
  );
});

d3.selectAll("input[name=dataset]").on("change", (e) => {
  let [xMin, xMax] = chromosomeView.getZoomRange();
  chromosomeView.update(
    cnvData[genomeView.getSelectedChromosome()],
    xMin,
    xMax
  );
  genomeView.update();
  if (resultsTable) {
    resultsTable.update();
  }
});
