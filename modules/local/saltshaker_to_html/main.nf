process SALTSHAKER_TO_HTML {
    tag "$meta.id"
    label "process_low"

    input:
    tuple val(meta), path(classify)

    output:
    tuple val(meta), path("*.html"), emit: classify_html

    script:
    """
    python3 << 'EOF'
    import re
    def saltshaker_txt_to_html(txt_file):
        with open(txt_file) as f:
            content = f.read()
        html_content = re.sub(r'\\t', '</td><td>', content)
        html_content = re.sub(r'\\n', '</td></tr><tr><td>', html_content)
        return html_content

    html = saltshaker_txt_to_html("${classify}")
    with open("${classify.baseName}.html", 'w') as f:
        f.write('<html><body><table border="1">')
        f.write(html)
        f.write('</table></body></html>')
    EOF
    """
}
