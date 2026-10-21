import base64
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from dash._callback_context import context_value
from dash._utils import AttributeDict

from autoqy_core.tools.smoother_gui import create_app
from autoqy_core.tools.spectral_recovery import SpectralRecovery, encode_session


def callback(app, part):
    return next(v['callback'].__wrapped__ for k,v in app.callback_map.items() if part in k)


class SpectralDialogTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.app = create_app(self.root / 'autosave.aqys')
        self.client = self.app.server.test_client()
        self.load = callback(self.app, 'load-message.children')
        self.save = callback(self.app, 'epsilon-save-message.children')
        self.reference = callback(self.app, 'epsilon-reference-upload.contents')

    def open_files(self, paths, data=None, initial=None):
        token = context_value.set(AttributeDict(triggered_inputs=[{'prop_id':'open-local-spectra.n_clicks'}], updated_props={}))
        try:
            with patch('autoqy_core.tools.smoother_gui._choose_files', return_value=paths) as dialog:
                result = self.load(None, 0, 1, [], [], [], None, data, None, [], [], [], [], [], [])
                dialog.assert_called_once_with(initial)
                return result
        finally:
            context_value.reset(token)

    def test_loading_reference_csv_image_and_session_share_latest_folder(self):
        dirs = []
        for name in ('spectra','csv','reference','image','session'):
            p=self.root/name
            p.mkdir()
            dirs.append(p)
        spectrum=dirs[0]/'input.csv'
        spectrum.write_text('Wavelength,value\n400,.2\n450,.4\n500,.1\n')
        data=self.open_files([str(spectrum)], initial=None)[0]
        destination=dirs[1]/'export.csv'
        with patch('autoqy_core.tools.smoother_gui._choose_save_file', return_value=str(destination)) as dialog:
            self.assertIn('Saved',self.save(1,None,data,[],source='old-folder')[0])
            self.assertEqual(dialog.call_args.args[1],str(dirs[0]))
        reference=dirs[2]/'epsilon.csv'
        reference.write_text('Wavelength_nm,Epsilon_mean_M-1_cm-1,Epsilon_SD_M-1_cm-1\n400,10,1\n450,20,1\n500,5,1\n')
        with patch('autoqy_core.tools.smoother_gui._choose_files',return_value=[str(reference)]) as dialog:
            content,filename,_=self.reference(1)
            dialog.assert_called_once_with(str(dirs[1]))
            self.assertEqual(filename,reference.name)
            self.assertEqual(base64.b64decode(content.split(',')[1]),reference.read_bytes())
        image=dirs[3]/'plot.svg'
        with patch('autoqy_core.tools.smoother_gui._choose_save_file',return_value=str(image)) as dialog:
            response=self.client.post('/_autoqy_spectral_image?filename=plot.svg',data=b'<svg xmlns="http://www.w3.org/2000/svg"/>',content_type='image/svg+xml')
            self.assertEqual(response.status_code,200,response.json)
            dialog.assert_called_once_with('plot.svg',str(dirs[2]))
        self.assertTrue(image.is_file())
        # Even opening a portable session must keep the folder just chosen, not
        # the original folder stored inside that session on another computer.
        session=dirs[4]/'saved.aqys'
        session.write_bytes(encode_session({'version':1,'dataset':data,'source':'unrelated-old-folder'}))
        restored=self.open_files([str(session)],data,initial=str(dirs[3]))
        self.assertIn('Restoring saved session',restored[1])
        self.open_files([],data,initial=str(dirs[4]))
        with patch('autoqy_core.tools.smoother_gui._choose_save_file',return_value=None) as dialog:
            self.assertIn('cancelled', self.save(2,None,data,[],source='old-folder')[0])
            self.assertEqual(dialog.call_args.args[1],str(dirs[4]))

    def test_image_cancel_and_invalid_payload_do_not_save_or_change_folder(self):
        with patch('autoqy_core.tools.smoother_gui._choose_save_file',return_value=None) as dialog:
            response=self.client.post('/_autoqy_spectral_image?filename=x.svg',data=b'<svg/>')
            self.assertEqual(response.json['message'],'Save cancelled.')
            dialog.assert_called_once_with('x.svg',None)
        with patch('autoqy_core.tools.smoother_gui._choose_save_file') as dialog:
            for name,payload in (('bad.png',b'not png'),('bad.txt',b'anything'),('../bad.svg',b'<svg/>')):
                response=self.client.post('/_autoqy_spectral_image?filename='+name,data=payload)
                self.assertEqual(response.status_code,400)
            dialog.assert_not_called()

    def test_image_folder_is_available_to_the_next_csv_save(self):
        path=self.root/'export.png'
        with patch('autoqy_core.tools.smoother_gui._choose_save_file',return_value=str(path)):
            response=self.client.post('/_autoqy_spectral_image?filename=export.png',data=b'\x89PNG\r\n\x1a\nexample')
            self.assertEqual(response.status_code,200)
        self.open_files([],initial=str(self.root))

    def test_last_dialog_folder_can_be_recovered(self):
        SpectralRecovery(self.root/'recover.aqys').write({'version':1,'source':'old-folder','dialog_folder':str(self.root)})
        app=create_app(self.root/'recover.aqys')
        load=callback(app,'load-message.children')
        token=context_value.set(AttributeDict(triggered_inputs=[{'prop_id':'open-local-spectra.n_clicks'}]))
        try:
            with patch('autoqy_core.tools.smoother_gui._choose_files',return_value=[]) as dialog:
                load(None,0,1,[],[],[],None,None,None,[],[],[],[],[],[])
                dialog.assert_called_once_with(str(self.root))
        finally:
            context_value.reset(token)

    def test_manual_session_save_uses_latest_folder_and_updates_it(self):
        from autoqy_core.smoother import SpectralDataset
        from autoqy_core.tools.smoother_gui import _pack
        from autoqy_core.tools.spectral_recovery import decode_session
        import numpy as np

        def components(node):
            if isinstance(node, (list, tuple)):
                for child in node:
                    yield from components(child)
            elif node is not None:
                yield node
                yield from components(getattr(node, "children", None))
        data = _pack(SpectralDataset(np.array([400.,450.,500.]),np.array([0.]),
                                    np.array([[.1],[.2],[.3]]),"csv"),["sample"],["sample.csv"])
        entry = next(v for k,v in self.app.callback_map.items() if 'session-save-message.children' in k)
        by_id = {getattr(c,"id",None): c for c in components(self.app.layout)
                 if isinstance(getattr(c,"id",None),str)}
        settings = [getattr(by_id[item["id"]], "value", None) for item in entry["state"][4:-6]]
        saved = self.root/'chosen'/'session.aqys'
        saved.parent.mkdir()
        with patch('autoqy_core.tools.smoother_gui._choose_save_file',return_value=str(saved)) as dialog:
            message = entry['callback'].__wrapped__(1,data,None,str(self.root),None,*settings,
                [1e-4],[1],[['on']],['sample'],
                [{'type':'legend-spectrum','index':0}],[{'type':'legend-name','index':0}])
            self.assertIn('Saved session',message)
            self.assertEqual(dialog.call_args.args[1],str(self.root))
        self.assertEqual(decode_session(saved.read_bytes())['dialog_folder'],str(saved.parent))
        self.open_files([],initial=str(saved.parent))

    def test_autosave_retains_the_shared_folder(self):
        image=self.root/'plot.svg'
        with patch('autoqy_core.tools.smoother_gui._choose_save_file',return_value=str(image)):
            self.client.post('/_autoqy_spectral_image?filename=plot.svg',data=b'<svg/>')
        response=self.client.post('/_autoqy_spectral_autosave',json={'version':1,'dataset':None})
        self.assertEqual(response.status_code,204)
        self.assertEqual(SpectralRecovery(self.root/'autosave.aqys').read()['dialog_folder'],str(self.root))
