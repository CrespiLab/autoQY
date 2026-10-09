/* Build spreadsheet clipboard text locally, including virtualized off-screen rows. */
window.autoqySheetClipboard = function (columns, rows, cells, copyAll, selectedRows, selectedColumns) {
    columns = columns || [];
    rows = rows || [];
    if (!columns.length || !rows.length) return "";
    const checkedRows = new Set((selectedRows || []).filter(row =>
        Number.isInteger(row) && row >= 0 && row < rows.length));
    const checkedColumns = new Set((selectedColumns || []).filter(id =>
        columns.some(column => column.id === id)));
    const quote = value => {
        if (value === null || value === undefined) return "";
        if (typeof value === "number" && !Number.isFinite(value)) return "";
        const text = String(value);
        return /[\t\r\n"]/.test(text) ? '"' + text.replace(/"/g, '""') + '"' : text;
    };
    const header = column => quote(Array.isArray(column.name) ? column.name.join(" · ") : column.name);
    if (!copyAll && (checkedRows.size || checkedColumns.size)) {
        const chosenColumns = columns.filter(column => !checkedColumns.size || checkedColumns.has(column.id));
        const lines = [chosenColumns.map(header).join("\t")];
        rows.forEach((row, index) => {
            if (!checkedRows.size || checkedRows.has(index)) {
                lines.push(chosenColumns.map(column => quote(row[column.id])).join("\t"));
            }
        });
        return lines.join("\r\n");
    }
    const selected = (cells || []).filter(cell =>
        Number.isInteger(cell.row) && cell.row >= 0 && cell.row < rows.length &&
        Number.isInteger(cell.column) && cell.column >= 0 && cell.column < columns.length &&
        (!cell.column_id || cell.column_id === columns[cell.column].id));
    if (!copyAll && !selected.length) return "";
    let firstRow = 0, lastRow = rows.length - 1;
    let firstColumn = 0, lastColumn = columns.length - 1;
    const selectedKeys = new Set();
    if (!copyAll) {
        firstRow = rows.length;
        lastRow = -1;
        firstColumn = columns.length;
        lastColumn = -1;
        for (const cell of selected) {
            firstRow = Math.min(firstRow, cell.row);
            lastRow = Math.max(lastRow, cell.row);
            firstColumn = Math.min(firstColumn, cell.column);
            lastColumn = Math.max(lastColumn, cell.column);
            selectedKeys.add(cell.row + ":" + cell.column);
        }
    }
    const chosenColumns = columns.slice(firstColumn, lastColumn + 1);
    const lines = [chosenColumns.map(header).join("\t")];
    for (let row = firstRow; row <= lastRow; row++) {
        lines.push(chosenColumns.map((column, offset) =>
            copyAll || selectedKeys.has(row + ":" + (firstColumn + offset))
                ? quote(rows[row][column.id]) : "").join("\t"));
    }
    return lines.join("\r\n");
};
