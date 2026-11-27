#!/usr/bin/env python3
import re
import typer
from pathlib import Path
from typing import List, Optional


CSS = """.tabcontent { display: none; }
.tabcontent.active { display: block; }
.tablinks { padding: 12px 20px; cursor: pointer; background-color: #f1f1f1; border: 1px solid #ccc; }
.tablinks.active { background-color: #4CAF50; color: white; }
.tablinks:hover { background-color: #ddd; }
"""


JS = """
function openTab(evt, tabName) {
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].classList.remove("active");
    }
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].classList.remove("active");
    }
    document.getElementById(tabName).classList.add("active");
    evt.currentTarget.classList.add("active");
}
document.getElementsByClassName("tablinks")[0].click();
"""

def sort_inputs(inputs, sample_ids):
    zipped = list(zip(inputs, sample_ids))
    sorted_zipped = sorted(zipped, key=lambda x: x[0].name)
    sorted_inputs = [x[0] for x in sorted_zipped]
    sorted_sample_ids = [x[1] for x in sorted_zipped]
    return sorted_inputs, sorted_sample_ids


def txt_to_html(txt_file):
    with open(txt_file) as tf:
        content = tf.read()
    html_content = re.sub(r'\\n', '<br>', content)
    return html_content


def create_tab_button(sample_id):
    return f'''<button class="tablinks" onclick="openTab(event, '{sample_id}')">{sample_id}</button>\n'''


def create_tab_content(sample_id, txt_file):
    html_content = txt_to_html(txt_file)
    return f'''<div id="{sample_id}" class="tabcontent">
\t<h3>{sample_id}</h3>
\t<pre style="padding: 15px; border-radius: 5px; overflow-x: auto;">{html_content}</pre>
</div>
'''


app = typer.Typer()

@app.command()
def main(
    input: List[Path] = typer.Option(
        ...,
        "--input",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Path to input .txt file (can be multiple)"
    ),
    sample: List[str] = typer.Option(
        ...,
        "--sample",
        help="Sample ID(s) corresponding to the input .txt file(s)"
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        help="Path to output .html file"
    )
):
    if len(input) != len(sample):
        raise typer.BadParameter(
            "--input and --sample must have the same number of values"
        )
    input, sample = sort_inputs(input, sample)
    tab_buttons = ''.join(create_tab_button(sid) for sid in sample)
    tab_contents = ''.join(create_tab_content(sid, inp) for inp, sid in zip(input, sample))

    html = f"""<html>
<head>
<style>
{CSS}
</style>
</head>
<body>
<div class="tab">
{tab_buttons}
</div>
{tab_contents}
<script>
{JS}
</script>
</body>
</html>"""

    with open(output, 'w') as f:
        f.write(html)

if __name__ == "__main__":
    app()
