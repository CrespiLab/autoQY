/* Run with node tests/test_sheet_clipboard.js. */
const assert = require("node:assert/strict");
global.window = {};
require("../autoqy_core/assets/sheet_clipboard.js");
const copy = window.autoqySheetClipboard;
const columns = [{id: "x", name: "Wavelength"}, {id: "a", name: "Edited legend"}, {id: "b", name: "PSS"}];
const rows = Array.from({length: 305}, (_, i) => ({x: 400 + i, a: i * 2, b: i * 3}));
assert.equal(copy(columns, rows, [], false), "");
assert.equal(copy(columns, rows, [], true).split("\r\n").length, 306);
assert.equal(copy(columns, rows, [], false, [304, 0], []),
    "Wavelength\tEdited legend\tPSS\r\n400\t0\t0\r\n704\t608\t912");
assert.equal(copy(columns, rows, [], false, [], ["b", "a"]).split("\r\n")[0], "Edited legend\tPSS");
assert.equal(copy(columns, rows, [], false, [2], ["a"]), "Edited legend\r\n4");
assert.equal(copy(columns, rows, [{row: 1, column: 1}, {row: 2, column: 1}], false), "Edited legend\r\n2\r\n4");
assert.equal(copy(columns, rows, [{row: 1, column: 0}, {row: 2, column: 1}], false),
    "Wavelength\tEdited legend\r\n401\t\r\n\t4");
assert.equal(copy(columns, rows, [{row: 1, column: 1, column_id: "old"}], false), "");
assert.equal(copy([{id: "a", name: 'Title\t"quoted"'}], [{a: null}, {a: 1.23456789012345}], [], true),
    '"Title\t""quoted"""\r\n\r\n1.23456789012345');
assert.equal(copy([], [], [], true), "");
console.log("Sheet clipboard tests passed (all rows, columns, intersection, ranges, escaping, precision).");
